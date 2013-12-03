////////////////////////////////////////////////////////////////////////////
//
// Copyright 1993-2013 NVIDIA Corporation.  All rights reserved.
//
// Please refer to the NVIDIA end user license agreement (EULA) associated
// with this source code for terms and conditions that govern your use of
// this software. Any use, reproduction, disclosure, or distribution of
// this software and related documentation outside the terms of the EULA
// is strictly prohibited.
//
////////////////////////////////////////////////////////////////////////////

/* Template project which demonstrates the basics on how to setup a project
 * example application.
 * Host code.
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

// includes CUDA
#include <cuda_runtime.h>

// includes, project
#include <helper_cuda.h>
#include <helper_functions.h> // helper functions for SDK examples
#include <device_functions.h> // helper functions for SDK examples



#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/gather.h>
#include <thrust/iterator/counting_iterator.h>

#include "libhough.h"

#undef OLD
#undef OLDMAP

typedef struct {
  unsigned int ntheta,nrho;
  float thetamin,thetamax;
  float rmin,rmax;
} houghLimits;


cudaStream_t streams[64];
void createStreams(unsigned int nb)
{
  if (nb>0 && nb<65)
    {
      for (unsigned int i=0;i<nb;i++)
	checkCudaErrors(cudaStreamCreate(&streams[i]));
    }
  
}

void deleteStreams(unsigned int nb)
{
  if (nb>0 && nb<65)
    {
      for (unsigned int i=0;i<nb;i++)
	checkCudaErrors(cudaStreamDestroy(streams[i]));
    }
  
}


//declare constant memory
__constant__ float limits[6];
__global__ void
checkLayerKernel(unsigned int *d_layer,unsigned int* d_pattern)
{
  __shared__ unsigned int sdata[1];
  //  const unsigned int nstub=blockDim.x; 
  const unsigned int istub=threadIdx.x;
  sdata[0]=0;
  __syncthreads();
  sdata[0] |=(1<<d_layer[istub]);
  __syncthreads();
  d_pattern[0]=sdata[0];
}
__global__ void
cleanUIKernel(unsigned int *d_layer)
{
  d_layer[blockIdx.x]=0;
}
__global__ void
cleanUI1Kernel(unsigned int *d_hough,unsigned int *d_hough_layer,unsigned int *d_hough_map,houghLimits hl)
{
  const unsigned int ith=threadIdx.x;
  const unsigned int ir=blockIdx.x;
  //const unsigned int nbintheta=int(limits[4]);
  const unsigned nbinrho=hl.nrho;//int(limits[5]);
  d_hough_layer[ith*nbinrho+ir]=0;
  d_hough[ith*nbinrho+ir]=0;
#ifdef OLDMAP
  for (int i=0;i<GPU_MAX_STUB_BIN;i++)
    d_hough_map[ith*nbinrho*GPU_MAX_STUB_BIN+ir*GPU_MAX_STUB_BIN+i]=0;
#endif
}
__global__ void
cleanFKernel(float *d_layer)
{
  d_layer[blockIdx.x]=0;
}

__device__ void 
localRegression(float xp,float yp,float* d,float* res)
{

  double z2=d[1],z=d[0],zx=d[2],x=d[3],n=d[4];
  z/=n;
  z2/=n;
  zx/=n;
  x/=n;
  double s2z = z2-z*z;
  double szx = zx-z*x;
  
  double a = szx/s2z;
  double b=x -a*z;
  double phi=atan(a);
  double theta=phi+PI/2.;
  double r=b*sin(theta);
  double R=0.5/abs(r);
  double pt=0.3*3.8*R/100.;
  double xi=-a/2./b;
  //  double yi=1./2./b;
  if (phi<0) phi+=2*PI;
  if (xp>0 && yp>0 && phi>PI) phi-=PI;
  if (xp<0 && yp>0 && phi>PI) phi-=PI;
  if (xp<0 && yp<0 && phi<PI) phi+=PI;
  if (xp>0 && yp<0 && phi<PI) phi+=PI;
  


  res[0]=a;
  res[1]=b;
  res[2]=phi;
  res[3]=theta;
  res[4]=r;
  res[5]=R;
  res[6]=pt;
  res[7]=xi;
  res[8]=-log(abs(tan(atan(a)/2)));
  if (z<0) res[8]=-res[8];
  res[9]=n;

 
}
__global__ void
doRegressionKernel(float *d_x, float *d_y,unsigned int* d_layer,unsigned int mode,float* d_reg)
{
  // shared memory
  // the size is determined by the host application
  //const unsigned int nbintheta=blockDim.x;

  __shared__  double fdata[5*1024];

  // access thread id
  const unsigned int id = threadIdx.x;
    
  // write data to global memory
  for (int iw=0;iw<5*1024;iw++)
    fdata[iw]=0;
  __syncthreads();

  unsigned int l=d_layer[id];
  if ((mode==1 && (l==5||l==6||l==7)) || mode==0)
    //  if ((1>0))
    { 
      fdata[5*id+0]+=d_x[id];
      fdata[5*id+1]+=d_x[id]*d_x[id];
      fdata[5*id+2]+=d_x[id]*d_y[id];
      fdata[5*id+3]+=d_y[id];
      fdata[5*id+4]+=1.;
    }
  //sdata[tid*nrhoword]=12;
  __syncthreads();
  __threadfence();
  if (id==1)
    {
      double z2=0,z=0,zx=0,x=0,n=0;
      for (int iw=0;iw<blockDim.x;iw++)
	{
	      z+=fdata[5*iw+0];
	      z2+=fdata[5*iw+1];
	      zx+=fdata[5*iw+2];
	      x+=fdata[5*iw+3];
	      n+=fdata[5*iw+4];

	}
      z/=n;
      z2/=n;
      zx/=n;
      x/=n;
      double s2z = z2-z*z;
      double szx = zx-z*x;
      
      double a = szx/s2z;
      double b=x -a*z;
      double phi=atan(a);
      double theta=phi+PI/2.;
      double r=b*sin(theta);
      double R=0.5/abs(r);
      double pt=0.3*3.8*R/100.;
      double xi=-a/2./b;
      double yi=1./2./b;
      if (phi<0) phi+=2*PI;
      if (d_x[0]>0 && d_y[0]>0 && phi>PI) phi-=PI;
      if (d_x[0]<0 && d_y[0]>0 && phi>PI) phi-=PI;
      if (d_x[0]<0 && d_y[0]<0 && phi<PI) phi+=PI;
      if (d_x[0]>0 && d_y[0]<0 && phi<PI) phi+=PI;



      d_reg[0]=a;
      d_reg[1]=b;
      d_reg[2]=phi;
      d_reg[3]=theta;
      d_reg[4]=r;
      d_reg[5]=R;
      d_reg[6]=pt;
      d_reg[7]=xi;
      d_reg[8]=yi;
      d_reg[9]=n;
    }
  __syncthreads();
  __threadfence();

}


__global__ void
calculateHoughPointKernel(float *d_x, float *d_y,unsigned int nbinrho,unsigned int* d_images)
{
  // shared memory
  // the size is determined by the host application
  //const unsigned int nbintheta=blockDim.x;
  const unsigned int nrhoword=nbinrho/32+1;
  __shared__  unsigned int sdata[512*16];

  // access thread id
  const unsigned int tid = threadIdx.x;
    
  // write data to global memory
  for (int iw=0;iw<nrhoword;iw++)
    sdata[tid*nrhoword+iw]=0;
  __syncthreads();

  
  double theta=limits[0]+(tid+0.5)*(limits[1]-limits[0])/blockDim.x;
  double r=d_x[blockIdx.x]*cos(theta)+d_y[blockIdx.x]*sin(theta);
  int ir=int(floor((r-limits[2])/((limits[3]-limits[2])/nbinrho)));
  if (ir>=0 && ir<nbinrho) 
    {
      int iw=ir/32;
      int ib=ir%32;
      sdata[tid*nrhoword+iw] |=(1<<ib);
    }
  //sdata[tid*nrhoword]=12;
  __syncthreads();

 
  // write data to global memory
  for (int iw=0;iw<nrhoword;iw++)
    d_images[blockIdx.x*nrhoword*blockDim.x+tid*nrhoword+iw] = sdata[tid*nrhoword+iw];
  //d_images[blockIdx.x]=blockIdx.x;
}

__global__ void
computeHoughPointKernel(float *d_x, float *d_y,short* d_val,houghLimits hl)
{
  // shared memory
  // the size is determined by the host application
  //const unsigned int nbintheta=blockDim.x;


  const unsigned int is=blockIdx.x;
  const unsigned int ith=threadIdx.x;
  const unsigned int nbintheta=hl.ntheta;//int(limits[4]);
  const unsigned nbinrho=hl.nrho;//int(limits[5]);
  const float thmin=hl.thetamin;
  const float thmax=hl.thetamax;
  const float rmin=hl.rmin;
  const float rmax=hl.rmax;
  
  
  double theta=thmin+(ith+0.5)*(thmax-thmin)/blockDim.x;
  double r=d_x[is]*cos(theta)+d_y[is]*sin(theta);
  short ir=int(floor((r-rmin)/((rmax-rmin)/nbinrho)));
  if (ir>=0 && ir<nbinrho) 
    {
      d_val[is*nbintheta+ith]=ir;
    }
  else
    d_val[is*nbintheta+ith]=-1;
  //sdata[tid*nrhoword]=12;
  __syncthreads();
}

__global__ void
fillHoughKernel(short *d_val,unsigned int* d_layer,unsigned int* d_hough,unsigned int* d_hough_layer,unsigned int* d_hough_map,houghLimits hl)
{
  // shared memory
  // the size is determined by the host application
  //const unsigned int nbintheta=blockDim.x;


  const unsigned int is=blockIdx.x;
  const unsigned int ith=threadIdx.x;
  const unsigned int nbinrho=hl.nrho;//int(limits[5]);
  const unsigned int nbintheta=hl.ntheta;//int(limits[4]);
  short ir= d_val[is*nbintheta+ith];
  if (ir>=0)
    {
      atomicAdd(&d_hough[ith*nbinrho+ir],1);
      atomicOr(&d_hough_layer[ith*nbinrho+ir],(1<<d_layer[is]));
#ifdef OLDMAP
      int iwm=is/32;
      int ibm=is%32;
      atomicOr(&d_hough_map[ith*nbinrho*GPU_MAX_STUB_BIN+ir*GPU_MAX_STUB_BIN+iwm],(1<<ibm));
#endif
    }

  __syncthreads();
}




__device__ unsigned int pointerIndex=0;
__global__ void
sumHoughPointKernel(unsigned int nstub,unsigned int* d_images,unsigned int* d_hough,unsigned int* d_hough_map,unsigned int* d_layer,unsigned int* d_hough_layer,unsigned int mode)
{
  const unsigned int ith=threadIdx.x;
  const unsigned int nth=blockDim.x;
  const unsigned int ir= blockIdx.x;
  const unsigned int nr= gridDim.x;
  double PT=1000.;
  if (mode == 0)
    {
      float r=limits[2]+(ir+0.5)*(limits[3]-limits[2])/nr;
      
      PT=1./2./fabs(r)*0.3*3.8/100.;
      d_hough[ith*nr+ir]=0;
      d_hough_layer[ith*nr+ir]=0;
#ifdef OLDMAP
      for (int iwm=0;iwm<GPU_MAX_STUB_BIN;iwm++)
	d_hough_map[ith*nr*GPU_MAX_STUB_BIN+ir*GPU_MAX_STUB_BIN+iwm]=0;
#endif
    }
  if (PT>1.4 ) {

  const unsigned int nrhoword=nr/32+1;
  const unsigned int iw=ir/32;
  const unsigned int ib=ir%32;
  d_hough[ith*nr+ir]=0;
  d_hough_layer[ith*nr+ir]=0;
  for (int iwm=0;iwm<GPU_MAX_STUB_BIN;iwm++)
    d_hough_map[ith*nr*GPU_MAX_STUB_BIN+ir*GPU_MAX_STUB_BIN+iwm]=0;
  for (int is=0;is<nstub;is++)
    {
      if (d_images[is*nrhoword*nth+ith*nrhoword+iw] & (1<<ib)) 
	{
	  d_hough[ith*nr+ir]+=1;
	  d_hough_layer[ith*nr+ir]|=(1<<d_layer[is]);
#ifdef OLDMAP
	  int iwm=is/32;
	  int ibm=is%32;
	  d_hough_map[ith*nr*GPU_MAX_STUB_BIN+ir*GPU_MAX_STUB_BIN+iwm]|= (1<<ibm);
#endif
	}
    }
  pointerIndex=0;
  }
}
__device__ unsigned int d_nbins_h=0;
__device__ float d_sum_h=0;
__device__ float d_sum2_h=0;

__global__ void
summaryHoughKernel(unsigned int* d_hough,float* d_param)
{
  unsigned int ith=threadIdx.x;
  unsigned int ir= blockIdx.x;
   if (ith ==1 && ir==1)
     {
     d_nbins_h=0;
     d_sum_h=0;
     d_sum2_h=0;
     }
   __syncthreads();
   __threadfence(); 
  if (d_hough[ith*gridDim.x+ir]>0)
    {
      atomicAdd(&d_nbins_h,1);
      atomicAdd(&d_sum_h,1.*d_hough[ith*gridDim.x+ir]);
      atomicAdd(&d_sum2_h,1.*d_hough[ith*gridDim.x+ir]*d_hough[ith*gridDim.x+ir]);
    }
  //if (ith==10 && ir==10) d_cand[0]=ith*gridDim.x+ir;
  __syncthreads();
  __threadfence();
  //if (ith==1 && ir==1)
  d_param[0]=d_nbins_h;
  d_param[1]=d_sum_h;
  d_param[2]=d_sum2_h;

}


__global__ void
ListHoughPointKernel(unsigned int* d_hough,unsigned int* d_hough_layer,unsigned int min_val,unsigned int min_layer,unsigned int* d_cand,houghLimits hl)
{
  const unsigned int nbinrho=hl.nrho;//int(limits[5]);
  const unsigned int nbintheta=hl.ntheta;//int(limits[4]);

  const unsigned int ith=threadIdx.x;
  const unsigned int ir= blockIdx.x;
   if (ith ==1 && ir==1)
     pointerIndex=0;
   __syncthreads();
   __threadfence(); 
  if (d_hough[ith*nbinrho+ir]>=min_val)
    {
      /*
      bool ismax=true;
      for (int ic=-1;ic<=1;ic++)
	for (int jc=-1;jc<=1;jc++)
	  {
	    if (ith==0 && ic<0) continue;
	    if (ith+ic>=nbintheta) continue;
	    if (ir==0 && jc<0) continue;
	    if (ir+jc>=nbinrho) continue;
	    if (d_hough[(ith+ic)*nbinrho+(ir+jc)]>d_hough[ith*nbinrho+ir])
	      {ismax=false;break;}
	  }
      */
      bool nmax=false;
      float val=d_hough[ith*nbinrho+ir]*1.;
      if (ith>0 && ir>0 && ith<=(nbintheta-1) && ir<(nbinrho-1))
	{
	  if ((val-d_hough[(ith-1)*nbinrho+ir])<0) nmax=true;
	  if ((val-d_hough[(ith+1)*nbinrho+ir])<0) nmax=true;
	  if((val-d_hough[(ith-1)*nbinrho+(ir-1)])<0) nmax=true;
	  if ((val-d_hough[(ith)*nbinrho+ir-1])<0) nmax=true;
	  if((val-d_hough[(ith+1)*nbinrho+ir-1])<0) nmax=true;
	  if((val-d_hough[(ith-1)*nbinrho+(ir+1)])<0) nmax=true;
	  if((val-d_hough[(ith)*nbinrho+ir+1])<0) nmax=true;
	  if((val-d_hough[(ith+1)*nbinrho+ir+1])<0) nmax=true;
	}
      if (!nmax)
	{
	  unsigned int pattern=d_hough_layer[ith*nbinrho+ir];

	  if (ith>0 && ir>0 && ith<=(nbintheta-1) && ir<(nbinrho-1))
	    {
	      pattern |=d_hough_layer[(ith-1)*nbinrho+ir];
	      pattern |=d_hough_layer[(ith+1)*nbinrho+ir];
	      pattern |=d_hough_layer[(ith-1)*nbinrho+ir-1];
	      pattern |=d_hough_layer[ith*nbinrho+ir-1];
	      pattern |=d_hough_layer[(ith+1)*nbinrho+ir-1];
	      pattern |=d_hough_layer[(ith-1)*nbinrho+ir+1];
	      pattern |=d_hough_layer[ith*nbinrho+ir+1];
	      pattern |=d_hough_layer[(ith+1)*nbinrho+ir+1];
	    }
	  pattern=d_hough_layer[ith*nbinrho+ir]; //@essai
	  unsigned int np=0;
	  bool l[24];
	  for (int ip=1;ip<=24;ip++)
	    {
	      l[ip]=((pattern &(1<<ip))!=0);
	      if (l[ip]) np++;
	    }
	  bool bar56=(l[5]&&l[6])||(l[5]&&l[7])||(l[6]&&l[7]);
	  // bar56=true;
	  //np=10;
	  if (np>=min_layer && d_hough[ith*nbinrho+ir]>=min_val && bar56)
	    {
	      unsigned int id=atomicInc(&pointerIndex,GPU_MAX_CAND);
	      //d_cand[0]+=1;
	      if (id<GPU_MAX_CAND-1)
		d_cand[id+1]=(ith & 0x3FF)+((ir&0x3FF)<<10)+((d_hough[ith*nbinrho+ir]&0x3FF)<<20);
	    }
	}
    }
  //if (ith==10 && ir==10) d_cand[0]=ith*gridDim.x+ir;
  __syncthreads();
  __threadfence();
  //if (ith==1 && ir==1)
  d_cand[0]=pointerIndex;

}

__global__ void
conformalPositionKernel(float* d_xo,float* d_yo,float* d_ro)
{
  const unsigned int ib=blockIdx.x;

  double r2=d_xo[ib]*d_xo[ib]+d_yo[ib]*d_yo[ib];
  double r=sqrt(r2);
  double x= d_xo[ib]/r2;
  double y= d_yo[ib]/r2;
  d_xo[ib]=x;
  d_yo[ib]=y;
  d_ro[ib]=r;
  __syncthreads();
  __threadfence();

}
__global__ void
copyNeighbourKernel(unsigned int ith,unsigned int ir,unsigned int ntheta,unsigned int nrho,unsigned int* d_hough_map,unsigned int* d_temp)
{
  if (threadIdx.x==1)
    {

  int ithmin=ith-1;
  int ithmax=ith+1;
  int irmin=ir-1;
  int irmax=ir+1;
  if (ith==0) ithmin=0;
  if (ir==0) irmin=0;
  if (ith==(ntheta-1)) ithmax=(ntheta-1);
  if (ir==(nrho-1)) ithmax=(nrho-1);
  d_temp[0]=0;
  for (int i=0;i<GPU_MAX_STUB_BIN;i++)
    {
      d_temp[i+1]=0;
    }
  for (int ii=ithmin;ii<=ithmax;ii++)
    for (int jj=irmin;jj<=irmax;jj++)
      {
	for (int iw=0;iw<GPU_MAX_STUB_BIN;iw++)
	  {

	  d_temp[iw+1]|=d_hough_map[ii*nrho*GPU_MAX_STUB_BIN+jj*GPU_MAX_STUB_BIN+iw];
	  }
      }
  for (int i=0;i<GPU_MAX_STUB_BIN;i++)
    for (int ib=0;ib<32;ib++)
      if (d_temp[i+1]&(1<<ib)) d_temp[0]+=1;
    
      
  __syncthreads();
    }
}

__global__ void
copyFromValKernel(unsigned int ith,unsigned int ir,unsigned int nbintheta,short* d_val,float* d_xi,float* d_yi,unsigned int* di_layer,float* d_ri,float* d_zi,float* d_xo,float* d_yo,unsigned int* do_layer,float* d_ro,float* d_zo,float* d_reg,bool regression,unsigned int* d_temp)
{
  const unsigned int ib=threadIdx.x;

  //__threadfence();
  if (d_val[ib*nbintheta+ith]==ir )
    {
      int iwm=ib/32;
      int ibm=ib%32;
      if (!(d_temp[iwm] & (1<<ibm)))
	{
	  d_temp[iwm]|=(1<<ibm); // no problem done bin/bin so one stub cannot be set in //
      float fid=atomicAdd(&d_reg[20],1.);
      unsigned int id=int(fid);
      //d_cand[0]+=1;
      if (id<GPU_MAX_STUB)
	{
	  float x=d_xi[ib],y=d_yi[ib],r=d_ri[ib],z=d_zi[ib];
	  unsigned int l=di_layer[ib]; 
	  d_xo[id]=x;
	  d_yo[id]=y;
	  d_ro[id]=r;
	  d_zo[id]=z;
	  do_layer[id]=l;
	  if (regression)
	    {
	      atomicAdd(&d_reg[50],x);
	      atomicAdd(&d_reg[51],x*x);
	      atomicAdd(&d_reg[52],x*y);
	      atomicAdd(&d_reg[53],y);
	      atomicAdd(&d_reg[54],1.);
	      if ((l==5) || (l==6) || (l==7))
		{
		  atomicAdd(&d_reg[55],z);
		  atomicAdd(&d_reg[56],z*z);
		  atomicAdd(&d_reg[57],z*r);
		  atomicAdd(&d_reg[58],r);
		  atomicAdd(&d_reg[59],1.);
		}
	    }
	}
	}
    }
  //if (ith==10 && ir==10) d_cand[0]=ith*gridDim.x+ir;
  __syncthreads();
  //__threadfence();
  //if (ith==1 && ir==1)
  //  d_cand[0]=pointerIndex;
  if (ib==1 && regression)
    {
      localRegression(d_xo[0],d_yo[0],&d_reg[50],&d_reg[60]);
      localRegression(d_xo[0],d_yo[0],&d_reg[55],&d_reg[70]);
    }
}


__global__ void
clearFloatKernel(float* d_float)
{
  d_float[threadIdx.x]=0;
  __syncthreads();
}
__global__ void
clearUIKernel(unsigned int* d_float)
{
  d_float[threadIdx.x]=0;
  __syncthreads();
}

__global__ void
copyPositionKernel(unsigned int* d_map,float* d_xi,float* d_yi,unsigned int* di_layer,float* d_ri,float* d_zi,float* d_xo,float* d_yo,unsigned int* do_layer,float* d_ro,float* d_zo,float* d_reg,bool regression)
{
  const unsigned int ib=threadIdx.x;

  if (ib==1)
    pointerIndex=0;
  for (int i=50;i<60;i++)
    d_reg[i]=0;
  __syncthreads();
  //__threadfence();
  int iwm=ib/32;
  int ibm=ib%32;
    if ((d_map[iwm]&(1<<ibm)) )
    {
      unsigned int id=atomicInc(&pointerIndex,512);
      //d_cand[0]+=1;
      if (id<512)
	{
	  float x=d_xi[ib],y=d_yi[ib],r=d_ri[ib],z=d_zi[ib];
	  unsigned int l=di_layer[ib]; 
	  d_xo[id]=x;
	  d_yo[id]=y;
	  d_ro[id]=r;
	  d_zo[id]=z;
	  do_layer[id]=l;
	  if (regression)
	    {
	      atomicAdd(&d_reg[50],x);
	      atomicAdd(&d_reg[51],x*x);
	      atomicAdd(&d_reg[52],x*y);
	      atomicAdd(&d_reg[53],y);
	      atomicAdd(&d_reg[54],1.);
	      if ((l==5) || (l==6) || (l==7))
		{
		  atomicAdd(&d_reg[55],z);
		  atomicAdd(&d_reg[56],z*z);
		  atomicAdd(&d_reg[57],z*r);
		  atomicAdd(&d_reg[58],r);
		  atomicAdd(&d_reg[59],1.);
		}
	    }
	}
    }
  //if (ith==10 && ir==10) d_cand[0]=ith*gridDim.x+ir;
  __syncthreads();
  //__threadfence();
  //if (ith==1 && ir==1)
  //  d_cand[0]=pointerIndex;
  if (ib==1 && regression)
    {
      localRegression(d_xo[0],d_yo[0],&d_reg[50],&d_reg[60]);
      localRegression(d_xo[0],d_yo[0],&d_reg[55],&d_reg[70]);
    }
}
void* h_malloc(unsigned int memSize)
{
  void* f;
  checkCudaErrors(cudaHostAlloc((void **)&f, memSize,cudaHostAllocWriteCombined));
  return f;
}
void h_free(void* f)
{
  checkCudaErrors(cudaFreeHost(f));
}


void createHough(houghParam* p)
{
   p->h_cand = (unsigned int*) malloc(GPU_MAX_CAND*sizeof(unsigned int));
   p->h_temp = (unsigned int*) malloc(512*sizeof(unsigned int));
   p->h_reg = (float*) malloc(GPU_MAX_REG*sizeof(float));
   p->h_val = (short*) malloc(GPU_MAX_STUB*GPU_MAX_THETA*sizeof(short));
   
   checkCudaErrors(cudaMalloc((void **) &p->d_reg,GPU_MAX_REG *sizeof(float)));
   checkCudaErrors(cudaMalloc((void **) &p->d_temp,512 *sizeof(unsigned int)));

   checkCudaErrors(cudaMalloc((void **) &p->d_val,GPU_MAX_STUB*GPU_MAX_THETA*sizeof(short)));
   checkCudaErrors(cudaMalloc((void **) &p->d_x,GPU_MAX_STUB *sizeof(float)));
   checkCudaErrors(cudaMalloc((void **) &p->d_y, GPU_MAX_STUB*sizeof(float)));
   checkCudaErrors(cudaMalloc((void **) &p->d_r, GPU_MAX_STUB*sizeof(float)));
   checkCudaErrors(cudaMalloc((void **) &p->d_z, GPU_MAX_STUB*sizeof(float)));
   checkCudaErrors(cudaMalloc((void **) &p->d_layer, GPU_MAX_STUB*sizeof(unsigned int)));
   checkCudaErrors(cudaMalloc((void **) &p->d_cand, GPU_MAX_CAND*sizeof(unsigned int)));
#ifdef OLD
   checkCudaErrors(cudaMalloc((void **) &p->d_images,GPU_MAX_STUB*GPU_MAX_THETA*GPU_MAX_RHO_WORD*sizeof(unsigned int)));
#endif
   checkCudaErrors(cudaMalloc((void **) &p->d_hough,GPU_MAX_THETA*GPU_MAX_RHO*sizeof(unsigned int) ));
   checkCudaErrors(cudaMalloc((void **) &p->d_hough_layer,GPU_MAX_THETA*GPU_MAX_RHO*sizeof(unsigned int) ));
#ifdef OLDMAP
   checkCudaErrors(cudaMalloc((void **) &p->d_hough_map,GPU_MAX_THETA*GPU_MAX_RHO*GPU_MAX_STUB_BIN*sizeof(unsigned int) ));
#endif

}
void clearHough(houghParam* p)
{
  //cleanUI1Kernel<<<p->nrho,p->ntheta>>>(p->d_hough,p->d_hough_layer,p->d_hough_map);
  //cleanUI1Kernel<<<GPU_MAX_RHO,GPU_MAX_THETA>>>(p->d_hough_layer);
}

void clean(houghParam* p)
{

  cleanFKernel<<<GPU_MAX_STUB,1>>>(p->d_x);
  cleanFKernel<<<GPU_MAX_STUB,1>>>(p->d_y);
  cleanFKernel<<<GPU_MAX_STUB,1>>>(p->d_z);
  cleanFKernel<<<GPU_MAX_STUB,1>>>(p->d_r);

  cleanUIKernel<<<GPU_MAX_STUB,1>>>(p->d_layer);
  cleanUIKernel<<<GPU_MAX_CAND,1>>>(p->d_cand);
#ifdef OLD
  cleanUIKernel<<<GPU_MAX_STUB*GPU_MAX_THETA*GPU_MAX_RHO_WORD,1>>>(p->d_images);
#endif
  cleanUIKernel<<<GPU_MAX_THETA*GPU_MAX_RHO,1>>>(p->d_hough);
  cleanUIKernel<<<GPU_MAX_THETA*GPU_MAX_RHO,1>>>(p->d_hough_layer);
  cleanUIKernel<<<GPU_MAX_STUB_BIN*GPU_MAX_THETA*GPU_MAX_RHO,1>>>(p->d_hough_map);
  cleanUIKernel<<<512,1>>>(p->d_temp);
}

void initialiseHough(houghParam* p,int nstub,int ntheta,int nrho,float thetamin,float thetamax,float rmin,float rmax)
{
  p->nstub=nstub;
  p->ntheta=ntheta;
  p->nrho=nrho;
  p->thetamin=thetamin;
  p->thetamax=thetamax;
  p->thetabin=(thetamax-thetamin)/ntheta;
  p->rmin=rmin;
  p->rmax=rmax;
  p->rbin=(rmax-rmin)/nrho;


}
void deleteHough(houghParam* p)
{
  free(p->h_cand);
  free(p->h_reg);
  free(p->h_temp);
  checkCudaErrors(cudaFree(p->d_val));
  checkCudaErrors(cudaFree(p->d_x));
  checkCudaErrors(cudaFree(p->d_y));
  checkCudaErrors(cudaFree(p->d_r));
  checkCudaErrors(cudaFree(p->d_z));
  checkCudaErrors(cudaFree(p->d_layer));
#ifdef OLD
  checkCudaErrors(cudaFree(p->d_images));
#endif
  checkCudaErrors(cudaFree(p->d_hough));
#ifdef OLDMAP
  checkCudaErrors(cudaFree(p->d_hough_map));
#endif
  checkCudaErrors(cudaFree(p->d_hough_layer));
  checkCudaErrors(cudaFree(p->d_cand));
  checkCudaErrors(cudaFree(p->d_reg));
  checkCudaErrors(cudaFree(p->d_temp));

}

void fillPositionHough(houghParam* p,float* h_x,float* h_y,float* h_z)
{
 checkCudaErrors(cudaMemcpy(p->d_x, h_x,p->nstub*sizeof(float),
			     cudaMemcpyHostToDevice));
 checkCudaErrors(cudaMemcpy(p->d_y, h_y,p->nstub*sizeof(float),
			     cudaMemcpyHostToDevice));
 checkCudaErrors(cudaMemcpy(p->d_z, h_z,p->nstub*sizeof(float),
			     cudaMemcpyHostToDevice));
}

void fillLayerHough(houghParam* p,unsigned int* h_layer)
{
 checkCudaErrors(cudaMemcpy(p->d_layer, h_layer,p->nstub*sizeof(unsigned int),
			     cudaMemcpyHostToDevice));
}

void copyHoughImage(houghParam* p,unsigned int* h_hough)
{
  checkCudaErrors(cudaMemcpy(h_hough, p->d_hough,p->ntheta*p->nrho*sizeof(unsigned int),
			     cudaMemcpyDeviceToHost));
 
}
void copyHoughLayer(houghParam* p,unsigned int* h_hough)
{
  checkCudaErrors(cudaMemcpy(h_hough, p->d_hough_layer,p->ntheta*p->nrho*sizeof(unsigned int),
			     cudaMemcpyDeviceToHost));
 
}

void fillConformalHough(houghParam* p,float* h_x,float* h_y,float* h_z)
{
  checkCudaErrors(cudaMemcpy(p->d_x, h_x,p->nstub*sizeof(float),
			     cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(p->d_y, h_y,p->nstub*sizeof(float),
			     cudaMemcpyHostToDevice));
 checkCudaErrors(cudaMemcpy(p->d_z, h_z,p->nstub*sizeof(float),
			     cudaMemcpyHostToDevice));
  dim3  grid1(p->nstub, 1, 1);
  conformalPositionKernel<<< grid1,1 >>>(p->d_x,p->d_y,p->d_r);
  cudaDeviceSynchronize();


}
void dump(houghParam* p)
{
  float x[1024];
  float y[1024];
  float z[1024];
  float r[1024];
  unsigned int layer[1024];
  
  checkCudaErrors(cudaMemcpy(x,p->d_x,p->nstub*sizeof(float),
			     cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(y,p->d_y,p->nstub*sizeof(float),
			     cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(z,p->d_z,p->nstub*sizeof(float),
			     cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(r,p->d_r,p->nstub*sizeof(float),
			     cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(layer,p->d_layer,p->nstub*sizeof(float),
			     cudaMemcpyDeviceToHost));

  for (int i=0;i<p->nstub;i++)
    printf("\t %d: (%f,%f,%f) r %f Layer %d \n",i,x[i],y[i],z[i],r[i],layer[i]);
}

void doRegression(houghParam* p,unsigned int mode )
{
  if (mode==0)
    doRegressionKernel<<<1,p->nstub>>>(p->d_x,p->d_y,p->d_layer,mode,p->d_reg);
  else
    if (mode==1)
      doRegressionKernel<<<1,p->nstub>>>(p->d_z,p->d_r,p->d_layer,mode,p->d_reg);
  checkCudaErrors(cudaMemcpy(p->h_reg,p->d_reg,GPU_MAX_REG*sizeof(float),
			     cudaMemcpyDeviceToHost));
}

void synchronize()
{
     cudaDeviceSynchronize();
     
}
void copyPositionHough(houghParam* pi,int icand,houghParam* po,unsigned int mode,bool regression,int streamid)
{
  cudaStream_t stream=0;
  if (streamid>=0) stream=streams[streamid]; 
   int ith=icand&0X3FF;
   int ir=(icand>>10)&0x3FF;

#ifdef OLDMAP
   dim3  grid1(1, 1, 1);
   dim3  grid2(512, 1, 1);

   // Copy neighbour maps to d_temp

   copyNeighbourKernel<<<2,128>>>(ith,ir,pi->ntheta,pi->nrho,pi->d_hough_map,pi->d_temp);
   cudaDeviceSynchronize();
   
   checkCudaErrors(cudaMemcpy(pi->h_temp,pi->d_temp,512*sizeof(unsigned int),
			      cudaMemcpyDeviceToHost));

   if (pi->h_temp[0]<((icand>>20)&0x3FF) || (mode==1))
     {
       //printf(" (%d,%d) Number of hits %d %d \n",ith,ir,(icand>>20)&0x3FF,pi->h_temp[0]);
       copyPositionKernel<<< 1,512 >>>(&pi->d_hough_map[ith*pi->nrho*GPU_MAX_STUB_BIN+ir*GPU_MAX_STUB_BIN],pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression);
   getLastCudaError("Kernel execution failed");
       po->nstub=(icand>>20)&0x3FF;
     }
   else
     {
       //   for (int k=0;k<GPU_MAX_STUB_BIN+1;k++)
       // printf("%d->%d %x\n",k,pi->h_temp[k],pi->h_temp[k]);
   //copyPositionKernel<<< 1,512 >>>(&pi->d_hough_map[ith*pi->nrho*GPU_MAX_STUB_BIN+ir*GPU_MAX_STUB_BIN],pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z);
       copyPositionKernel<<< 1,512 >>>(&pi->d_temp[1],pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression);
   //copyPositionKernel<<< 1,512 >>>(&pi->d_temp[1],pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z);
   cudaDeviceSynchronize();

   getLastCudaError("Kernel execution failed");
   cudaDeviceSynchronize();

   po->nstub=(icand>>20)&0x3FF;
   po->nstub=pi->h_temp[0];


   


     }
   if (regression)
     checkCudaErrors(cudaMemcpy(po->h_reg,po->d_reg,GPU_MAX_REG*sizeof(float),cudaMemcpyDeviceToHost));
#else
   //printf("Stream %d %x \n",streamid,(unsigned long) stream);
   clearFloatKernel<<<1,GPU_MAX_REG,0,stream>>>(po->d_reg);
   clearUIKernel<<<1,512,0,stream>>>(po->d_temp);

   if (mode == 1)
     {
       copyFromValKernel<<<1,pi->nstub,0,stream>>>(ith,ir,pi->ntheta,pi->d_val,pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression,po->d_temp);
     }
   else
     {
       copyFromValKernel<<<1,pi->nstub,0,stream>>>(ith,ir,pi->ntheta,pi->d_val,pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression,po->d_temp);
       if (streamid<0) cudaDeviceSynchronize();


       if (ith>0)
	 copyFromValKernel<<<1,pi->nstub,0,stream>>>(ith-1,ir,pi->ntheta,pi->d_val,pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression,po->d_temp);
       if (streamid<0) cudaDeviceSynchronize();


       if (ith<pi->ntheta-1)
	 copyFromValKernel<<<1,pi->nstub,0,stream>>>(ith+1,ir,pi->ntheta,pi->d_val,pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression,po->d_temp);
       if (streamid<0) cudaDeviceSynchronize();


       if (ith>0 && ir>0)
	 copyFromValKernel<<<1,pi->nstub,0,stream>>>(ith-1,ir-1,pi->ntheta,pi->d_val,pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression,po->d_temp);
       if (streamid<0) cudaDeviceSynchronize();


       if (ir>0)
	 copyFromValKernel<<<1,pi->nstub,0,stream>>>(ith,ir-1,pi->ntheta,pi->d_val,pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression,po->d_temp);
       if (streamid<0) cudaDeviceSynchronize();


       if (ir>0 && ith<pi->ntheta-1)
	 copyFromValKernel<<<1,pi->nstub,0,stream>>>(ith+1,ir-1,pi->ntheta,pi->d_val,pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression,po->d_temp);
       if (streamid<0) cudaDeviceSynchronize();


       if (ith>0 && ir<pi->nrho-1)
	 copyFromValKernel<<<1,pi->nstub,0,stream>>>(ith-1,ir+1,pi->ntheta,pi->d_val,pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression,po->d_temp);
       if (streamid<0) cudaDeviceSynchronize();


       if (ir<pi->nrho-1)
	 copyFromValKernel<<<1,pi->nstub,0,stream>>>(ith,ir+1,pi->ntheta,pi->d_val,pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression,po->d_temp);
       if (streamid<0) cudaDeviceSynchronize();


       if (ir<pi->nrho-1 &&  ith<pi->ntheta-1)
	 copyFromValKernel<<<1,pi->nstub,0,stream>>>(ith+1,ir+1,pi->ntheta,pi->d_val,pi->d_x,pi->d_y,pi->d_layer,pi->d_r,pi->d_z,po->d_x,po->d_y,po->d_layer,po->d_r,po->d_z,po->d_reg,regression,po->d_temp);
       if (streamid<0) cudaDeviceSynchronize();


     }
   getLastCudaError("Kernel execution failed");
   if (streamid<0) cudaDeviceSynchronize();
   checkCudaErrors(cudaMemcpyAsync(po->h_reg,po->d_reg,GPU_MAX_REG*sizeof(float),cudaMemcpyDeviceToHost,stream));
   if (streamid<0) 
     {
       cudaDeviceSynchronize();
       po->nstub=int(po->h_reg[20]);
     }
   //dump(po);
   //getchar();

#endif
   /*
   for (int i=50;i<60;i++)
     printf("%f ",po->h_reg[i]);
   printf("\n");

   for (int i=60;i<70;i++)
     printf("%f ",po->h_reg[i]);
   printf("\n");
   for (int i=70;i<80;i++)
     printf("%f ",po->h_reg[i]);
   printf("\n");


    getchar();
   */
   /*
   unsigned int* hmap= (unsigned int*) malloc(GPU_MAX_STUB_BIN*sizeof(unsigned int));
   checkCudaErrors(cudaMemcpy(hmap,&pi->d_hough_map[ith*pi->nrho*GPU_MAX_STUB_BIN+ir*GPU_MAX_STUB_BIN],GPU_MAX_STUB_BIN*sizeof(unsigned int),
			     cudaMemcpyDeviceToHost));

   printf("%x %d\n",icand,po->nstub);
   for (int iw=0;iw<GPU_MAX_STUB_BIN;iw++)
     printf("%x ",hmap[iw]);
   printf("\n");

   float* hx=(float*) malloc(po->nstub*sizeof(float));
   float* hy=(float*) malloc(po->nstub*sizeof(float));
   checkCudaErrors(cudaMemcpy(hx,po->d_x,po->nstub*sizeof(float),cudaMemcpyDeviceToHost));
   checkCudaErrors(cudaMemcpy(hy,po->d_y,po->nstub*sizeof(float),cudaMemcpyDeviceToHost));

   for (int ip=0;ip<po->nstub;ip++)
     printf("%f %f \n",hx[ip],hy[ip]);
   getchar();
   free(hmap);
   free(hx);
   free(hy);
   */
}

void processHough(houghParam* p,unsigned int min_cut,unsigned int min_layer,unsigned int mode,int streamid)
{
  cudaStream_t stream=0;
  if (streamid>=0) stream=streams[streamid]; 
  /****
  float hlim[6];
  hlim[0]=p->thetamin;
  hlim[1]=p->thetamax;
  hlim[2]=p->rmin;
  hlim[3]=p->rmax;
  hlim[4]=p->ntheta;
  hlim[5]=p->nrho;
  cudaMemcpyToSymbol(limits,  hlim,   sizeof(float)*6  );
  ***/
  houghLimits hl;
  hl.thetamin=p->thetamin;
  hl.thetamax=p->thetamax;
  hl.rmin=p->rmin;
  hl.rmax=p->rmax;
  hl.ntheta=p->ntheta;
  hl.nrho=p->nrho;
  
  // setup execution parameters
  dim3  grid1(p->nstub, 1, 1);
  dim3  threads1(p->ntheta, 1, 1);
  dim3 grid2(p->nrho,1,1);
  // printf("%d %d %d === %d %x %d \n",p->nstub,p->ntheta,p->nrho,streamid,(unsigned long) stream,mode);
  //getchar();
  if (mode==0)
    {
    computeHoughPointKernel<<< grid1, threads1,0,stream>>>(p->d_x,p->d_y,p->d_val,hl);
    }
  else
    if (mode==1)
      computeHoughPointKernel<<< grid1,threads1,0,stream>>>(p->d_z,p->d_r,p->d_val,hl);
  getLastCudaError("Kernel execution failed");
  if (streamid<0)
    cudaDeviceSynchronize();
  cleanUI1Kernel<<<p->nrho,p->ntheta,0,stream>>>(p->d_hough,p->d_hough_layer,p->d_hough_map,hl);
  getLastCudaError("Kernel execution failed");
  if (streamid<0) cudaDeviceSynchronize();
  /*
    checkCudaErrors(cudaMemcpy(p->h_val,p->d_val,GPU_MAX_STUB*GPU_MAX_THETA*sizeof(short),
  cudaMemcpyDeviceToHost));
  for (int ii=0;ii<p->nstub;ii++)
    {
      for (int jj=0;jj<p->ntheta;jj++)
  	printf("%d ",p->h_val[ii*p->ntheta+jj]);
      printf("\n");
    }
  getchar();
  */
  fillHoughKernel<<< grid1, threads1,0,stream>>>(p->d_val,p->d_layer,p->d_hough,p->d_hough_layer,p->d_hough_map,hl);
  getLastCudaError("Kernel execution failed");


  p->max_val=0;
  unsigned int threshold=4;
  //int(floor(m+3*rms));
  if (min_cut!=0 ) threshold=min_cut;
  //if (threshold<int(floor(p->max_val*0.5))) threshold=int(floor(p->max_val*0.5));
  //threshold=int(floor(m+3*rms));
   //printf("Max val %d Threshold %d \n",p->max_val,threshold);
  ListHoughPointKernel<<< grid2,threads1,0,stream>>>(p->d_hough,p->d_hough_layer,threshold,min_layer,p->d_cand,hl);
  getLastCudaError("Kernel execution failed");
  if (streamid<0) cudaDeviceSynchronize();

  checkCudaErrors(cudaMemcpyAsync(p->h_cand, p->d_cand, GPU_MAX_CAND*sizeof(unsigned int) ,
				  cudaMemcpyDeviceToHost,stream));
  
  }



StopWatchInterface *theTimer=0;

void initialiseTimer()
{
  theTimer=0;
  sdkCreateTimer(&theTimer);
}

void startTimer()
{
   sdkResetTimer(&theTimer);
   sdkStartTimer(&theTimer);
}

float stopTimer()
{
  sdkStopTimer(&theTimer);
  float t=sdkGetTimerValue(&theTimer);
  //printf("Processing time: %f (ms)\n",t );
  return t;
}
void deleteTimer()
{
   sdkDeleteTimer(&theTimer);
}
#ifdef RUN_TEST
void
runTest(int argc, char **argv)
{
 

  printf("%s Starting...\n\n", argv[0]);

  // use command-line specified CUDA device, otherwise use device with highest Gflops/s
  int devID = findCudaDevice(argc, (const char **)argv);

 


 
  StopWatchInterface *timer1 = 0;
  sdkCreateTimer(&timer1);
 

  
  unsigned int nstub=250;
  unsigned int ntheta=128;
  unsigned int nrho=256;
  //  unsigned int nrhow=nrho/32+1;
  
  // allocate host memory
  float *h_x = (float *) malloc(nstub*sizeof(float));
  float *h_y = (float *) malloc(nstub*sizeof(float));
  unsigned int* h_cand = (unsigned int*) malloc(GPU_MAX_CAND*sizeof(unsigned int));
  unsigned int* h_cand1 = (unsigned int*) malloc(GPU_MAX_CAND*sizeof(unsigned int));
  // initalize the memory
  for (unsigned int i = 0; i < nstub; ++i)
    {
      h_x[i] = (float) i*10./nstub;
      h_y[i] = (2*h_x[i]+1.);
    }
  float ymin=0,ymax=0;
  float totaltime=0;

  
  houghParam ph;
  createHough(&ph);


  houghParam phi;
  createHough(&phi);




  for (int jj=0;jj<1000;jj++)
    {
      nstub=rand()%(240-20)+20;
      for (unsigned int i = 0; i < nstub; ++i)
	{
	  
	  h_x[i] = rand()/((double)RAND_MAX ) * 10;
	  if (i%3!=0)
	    {
	      ymin=(1.56*h_x[i]+2.)*0.98;
	      ymax=(1.56*h_x[i]+2.)*1.02;
	    }
	  else
	    {
	      ymin=0;
	      ymax=21.;
	    }
	  
	  h_y[i] =  rand()/((double)RAND_MAX ) *(ymax-ymin)+ymin;
	}
      sdkResetTimer(&timer1);
      sdkStartTimer(&timer1);
      initialiseHough(&ph,nstub,ntheta,nrho,-PI/2,PI/2,-21.,21.);
      fillPositionHough(&ph,h_x,h_y);
      processHough(&ph,0);
      //doHough(nstub, h_x,h_y,ntheta,nrho,-PI/2.,PI/2,-21.,21.,h_cand);
      for (int ic=0;ic<ph.h_cand[0];ic++)
	{
	  int ith=ph.h_cand[ic+1]&0X3FF;
	  int ir=(ph.h_cand[ic+1]>>10)&0x3FF;

	  float tmi=-PI/2+(ith-1)*PI/ntheta;
	  float tma=-PI/2+(ith+1)*PI/ntheta;
	  float rmi=-21.+(ir-1)*42./nrho;
	  float rma=-21.+(ir+1)*42./nrho;
	  //printf("Look for bin %x  %d %d %f %f %f %f \n",ph.h_cand[ic+1],ith,ir,tmi,tma,rmi,rma);
	  //getchar();
	  initialiseHough(&phi,nstub,64,64,tmi,tma,rmi,rma);
	  copyPositionHough(&ph,ph.h_cand[ic+1],&phi);
	  processHough(&phi,0);
	  //doHough(nstub, h_x,h_y,64,64,tmi,tma,rmi,rma,h_cand1);

	  //getchar();
	   
	}
      cudaDeviceSynchronize();

      sdkStopTimer(&timer1);
      float curti=sdkGetTimerValue(&timer1);
      totaltime+=curti;
      if (jj%10==9)
      printf("%d %d Processing time: %f (ms) %f\n", jj,nstub,curti,totaltime*1E-3);
      
      continue;
    }
  sdkStopTimer(&timer1);
  printf("Processing time: %f (ms)\n", sdkGetTimerValue(&timer1));
  sdkDeleteTimer(&timer1);
  deleteHough(&phi);
  deleteHough(&ph);
  cudaDeviceReset();

}
#endif
