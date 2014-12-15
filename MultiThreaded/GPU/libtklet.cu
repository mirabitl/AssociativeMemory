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
#include "libtklet.h"

void createTkletEvent(cudatklevent* e)
{
  e->host_= (ctklevent*) malloc(sizeof(ctklevent));
 
  checkCudaErrors(cudaMalloc((void **) &e->device_,sizeof(ctklevent)));
}

void deleteTkletEvent(cudatklevent* e)
{
  free(e->host_);
  checkCudaErrors(cudaFree(e->device_));
}

void copyFromHost(cudatklevent* e)
{
  checkCudaErrors(cudaMemcpy(e->device_->flay_,e->host_->flay_,3*GPU_MAX_STUB*sizeof(float)+(32*128*sizeof(int)),cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(&e->device_->lay_,&e->host_->lay_,(GPU_MAX_STUB+5)*sizeof(unsigned int),cudaMemcpyHostToDevice));
}

void copyToHost(cudatklevent* e)
{
  checkCudaErrors(cudaMemcpy(&(e->host_->ntkl_),&(e->device_->ntkl_),sizeof(unsigned int)+GPU_MAX_CTKLET*sizeof(ctklet),
			     cudaMemcpyDeviceToHost));
}

__global__ void
fillConformalKernel(ctklevent* d)
{
  const unsigned int ib=blockIdx.x;

  double r2=d->x_[ib]*d->x_[ib]+d->y_[ib]*d->y_[ib];
  d->r_[ib]=__fsqrt_ru(r2);
  d->xp_[ib]=d->x_[ib]/r2;
  d->yp_[ib]=d->y_[ib]/r2;
  //__syncthreads();
  d->ntkl_=0;

}


void fillDevice(cudatklevent* e)
{
  copyFromHost(e);
  dim3  grid1(e->host_->nstub_, 1, 1);
  fillConformalKernel<<< grid1,1,0>>>(e->device_);

  cudaDeviceSynchronize();
}

__global__ void
combineLayerKernel(ctklevent* d)
{
  if (d->ntkl_==GPU_MAX_CTKLET) return;
  const unsigned int is=blockIdx.x;
  const unsigned int js=threadIdx.x;
 // Ask z info
  if ( ((d->lay_[is]>>16)&0x3)==0) return;
  if ( ((d->lay_[js]>>16)&0x3)==0) return;
  const unsigned int l1=(d->lay_[is]&0xFFFF);
  const unsigned int l2=(d->lay_[js]&0xFFFF);
  bool barok= ((l1==5 && l2==6) || (l1==6 && l2==7) || (l1==5 && l2==7));
  if (d->barrel_ && !barok) return;
  if (d->inter_ && !barok) return;
  if (d->endcap_ && d->sector_<8)
    {
      if (!((l1==5 && l2==6) || (l1==5 && l2==18) || (l1==5 && l2==19) || (l1==18 && l2==19)|| (l1==6 && l2==18) || (l1==6 && l2==19)  )) return;
    }
  if (d->endcap_ && d->sector_>=40)
    {
      if (!( (l1==5 && l2==6) || (l1==5 && l2==11) || (l1==5 && l2==12) || (l1==11 && l2==12)|| (l1==6 && l2==11) || (l1==6 && l2==12) )) return;
    }

  
 
  //
  double a =(d->yp_[js]-d->yp_[is])/(d->xp_[js]-d->xp_[is]);
  //if ((isel%4==0) && (a<-0.25 || a>1.55)) continue;
  double b =d->yp_[js]-a*d->xp_[js];
  double pt=5.7E-3*__fsqrt_ru((a*a+1)/b/b);
  if (fabs(pt)<1.8) return;
  double ar=(d->r_[js]-d->r_[is])/(d->z_[js]-d->z_[is]);
  double br=d->r_[js]-ar*d->z_[js];
  double zi=-br/ar;
  if (fabs(zi)>20) return;

  unsigned int id=atomicInc(&d->ntkl_,GPU_MAX_CTKLET);
  ctklet* tkl=&d->cand_[id];
	  //printf("%x %x %d %d \n",d->lay_[i],d->lay_[j],(d->lay_[i]>>16)&0x3,(d->lay_[j]>>16)&0x3);
  tkl->nhit_=2; 
  tkl->pattern_=0; 
  tkl->ok_=true;
  tkl->pattern_ |=(1<<l1);
  tkl->pattern_ |=(1<<l2);
  tkl->idx_[0]=is;
  tkl->idx_[1]=js;
    //__syncthreads();

}
__global__ void
combineLayerKernel(ctklevent* d,unsigned int l1,unsigned int l2)
{
  if (d->ntkl_==GPU_MAX_CTKLET) return;

  const unsigned int i1=blockIdx.x;
  if (i1>=128) return;
  const unsigned int j1=threadIdx.x;
  if (j1>=128) return;
  int is=d->flay_[l1*128+i1];
  if (is>GPU_MAX_STUB || is<0 ) return;
  int js=d->flay_[l2*128+j1]; 
  if (js>GPU_MAX_STUB || js<0) return;
  
  
  //if ((d->lay_[is]&0xFFFF)!=l1) return;
  //if ((d->lay_[js]&0xFFFF)!=l2) return;
 // Ask z info
  if ( ((d->lay_[is]>>16)&0x3)==0) return;
  if ( ((d->lay_[js]>>16)&0x3)==0) return;
  /* const unsigned int l1=(d->lay_[is]&0xFFFF); */
  /* const unsigned int l2=(d->lay_[js]&0xFFFF); */
  /* bool barok= ((l1==5 && l2==6) || (l1==6 && l2==7) || (l1==5 && l2==7)); */
  /* if (d->barrel_ && !barok) return; */
  /* if (d->inter_ && !barok) return; */
  /* if (d->endcap_ && d->sector_<8) */
  /*   { */
  /*     if (!( (l1==5 && l2==18) || (l1==5 && l2==19) || (l1==18 && l2==19))) return; */
  /*   } */
  /* if (d->endcap_ && d->sector_>=48) */
  /*   { */
  /*     if (!( (l1==5 && l2==11) || (l1==5 && l2==12) || (l1==11 && l2==12))) return; */
  /*   } */

  
 
  //
  double a =(d->yp_[js]-d->yp_[is])/(d->xp_[js]-d->xp_[is]);
  //if ((isel%4==0) && (a<-0.25 || a>1.55)) continue;
  double b =d->yp_[js]-a*d->xp_[js];
  double pt=5.7E-3*__fsqrt_ru((a*a+1)/b/b);
  if (fabs(pt)<1.8) return;
  double ar=(d->r_[js]-d->r_[is])/(d->z_[js]-d->z_[is]);
  double br=d->r_[js]-ar*d->z_[js];
  double zi=-br/ar;
  if (fabs(zi)>20) return;

  unsigned int id=atomicInc(&d->ntkl_,GPU_MAX_CTKLET);
  ctklet* tkl=&d->cand_[id];
	  //printf("%x %x %d %d \n",d->lay_[i],d->lay_[j],(d->lay_[i]>>16)&0x3,(d->lay_[j]>>16)&0x3);
  tkl->nhit_=2; 
  tkl->pattern_=0; 
  tkl->ok_=true;
  tkl->pattern_ |=(1<<l1);
  tkl->pattern_ |=(1<<l2);
  tkl->idx_[0]=is;
  tkl->idx_[1]=js;
    //__syncthreads();

}

void combineLayer(cudatklevent* e)
{
#undef OLDCOMB
#ifdef OLDCOMB
  dim3  grid1(e->host_->nstub_,1, 1);
  combineLayerKernel<<< grid1,e->host_->nstub_>>>(e->device_);
  cudaDeviceSynchronize();
#else
   
  if (e->host_->barrel_ || e->host_->inter_)
    { 
      dim3  grid1(128,1, 1);
      dim3  threads1(128,1, 1);
      combineLayerKernel<<<grid1,threads1>>>(e->device_,5,6);
      //cudaDeviceSynchronize();
      combineLayerKernel<<<128,128>>>(e->device_,5,7);
      //cudaDeviceSynchronize();
      combineLayerKernel<<<128,128>>>(e->device_,6,7);
      cudaDeviceSynchronize();
    }
  else
    {
      if (e->host_->sector_<8)
	{
      dim3  grid1(128,1, 1);
      dim3  threads1(128,1, 1);
      combineLayerKernel<<<grid1,threads1>>>(e->device_,5,6);
      combineLayerKernel<<<grid1,threads1>>>(e->device_,5,18);
      //cudaDeviceSynchronize();
      //combineLayerKernel<<<128,128>>>(e->device_,5,19);
      //cudaDeviceSynchronize();
      combineLayerKernel<<<grid1,threads1>>>(e->device_,6,18);
      //cudaDeviceSynchronize();
      //combineLayerKernel<<<128,128>>>(e->device_,6,19);

      //combineLayerKernel<<<128,128>>>(e->device_,18,19);
      cudaDeviceSynchronize();

	}
      else
	{
	        dim3  grid1(128,1, 1);
      dim3  threads1(128,1, 1);
      combineLayerKernel<<<grid1,threads1>>>(e->device_,5,6);
      combineLayerKernel<<<grid1,threads1>>>(e->device_,5,11);
      //cudaDeviceSynchronize();
      //combineLayerKernel<<<128,128>>>(e->device_,5,12);
      combineLayerKernel<<<grid1,threads1>>>(e->device_,6,11);
      //cudaDeviceSynchronize();
      // combineLayerKernel<<<128,128>>>(e->device_,6,12);

      //cudaDeviceSynchronize();
      // combineLayerKernel<<<128,128>>>(e->device_,11,12);
      cudaDeviceSynchronize();

	}
    }
#endif
}



__global__ void
addLayerKernel(ctklevent* e)
{
  const unsigned int it=blockIdx.x;
  const unsigned int is=threadIdx.x;
  ctklet* t=&(e->cand_[it]);
  if (!t->ok_) return;
  unsigned int l1=(e->lay_[is]&0xFFFF);

  if ((t->pattern_&(1<<l1))!=0) return; // layer already included

  //double r2=e->r_[is]*e->r_[is];
  //r2=1.;
  double distx=(t->ax_*e->xp_[is]+t->bx_-e->yp_[is])/__fsqrt_ru(1+t->ax_*t->ax_);
	  //	  if (fabs(distx)<1E-3)
	  //printf("is %d %f r2 %f %f %f %x pattern %x  \n",is,e->r_[is],distx,r2,t->ax_,e->lay_[is],t->pattern_);


	  //if (fabs(distx)>0.3) continue;
  
  double cut=3.5E-5;
  if (l1<=10) cut /=2;
  if (e->barrel_ && fabs(distx)>cut) return;
  if (e->inter_ && fabs(distx)>cut) return;
  if (e->endcap_ && fabs(distx)>0.7*cut) return;
  if (((e->lay_[is]>>16)&0x3)!=0)
    {
      double distr=(t->ar_*e->z_[is]+t->br_-e->r_[is])/__fsqrt_ru(1+t->ar_*t->ar_);
      
      if (!e->endcap_ && fabs(distr)>0.6) return;
	      //  if (e->barrel_ && fabs(distr)>0.45) continue;
	      // if (e->inter_ && fabs(distr)>0.2) continue;
      if (e->endcap_ && fabs(distr)>0.4) return;
    }
  unsigned int old=atomicOr(&t->pattern_,(1<<l1));
  if ((old&(1<<l1))!=0) return;
  unsigned int id=atomicInc(&t->nhit_,32);
  t->idx_[id]=is;

  //__syncthreads();
}    


void addLayer(cudatklevent* e)
{
  //copyToHost(e);
  checkCudaErrors(cudaMemcpy(&e->host_->ntkl_,&e->device_->ntkl_,sizeof(unsigned int),cudaMemcpyDeviceToHost));
  //printf("%d %d %d %d \n",e->host_->ntkl_,e->host_->endcap_,e->host_->barrel_,e->host_->inter_);
  // getchar();
  dim3  grid1(e->host_->nstub_,1, 1);
  dim3  grid2(e->host_->ntkl_,1, 1);
  
  addLayerKernel<<< grid2,grid1>>>(e->device_);
  cudaDeviceSynchronize();

}

__global__ void
computeTkletKernel(ctklevent* e)
{
  
  //const unsigned int it=blockIdx.x;
  const unsigned int it=threadIdx.x+1024*blockIdx.x;
  if ((it+1)>e->ntkl_) return;
  ctklet* t=&(e->cand_[it]);
  if (!t->ok_) return;
  t->pattern_=0;
  t->ok_=true;
  memset(&t->sumx_,0,20*sizeof(double));
  //printf("clearing %d %d \n",it,t->nhit_);
  int idx;
  for (int ih=0;ih<t->nhit_;ih++)
    {
      idx=t->idx_[ih];
      int layerh=(e->lay_[idx]&0xFFFF);
      t->pattern_ |=(1<<layerh);
      t->sumx_+=e->xp_[idx];
      t->sumx2_+=e->xp_[idx]*e->xp_[idx];
      t->sumxy_+=e->xp_[idx]*e->yp_[idx];
      t->sumy_+=e->yp_[idx];
      t->nxy_+=1.;
      if (((e->lay_[idx]>>16)&0x3)==0) continue;
      t->sumz_+=e->z_[idx];
      t->sumz2_+=e->z_[idx]*e->z_[idx];
      t->sumzr_+=e->z_[idx]*e->r_[idx];
      t->sumr_+=e->r_[idx];
      t->nzr_+=1.;

    }
  if (t->nzr_<2 || t->nxy_<2) 
    {
      t->ok_=false;
      return;
    }
  double s2z = t->sumx2_/t->nxy_-(t->sumx_/t->nxy_)*(t->sumx_/t->nxy_);
  double szx = t->sumxy_/t->nxy_-(t->sumx_/t->nxy_)*(t->sumy_/t->nxy_);
  
  t->ax_ = szx/s2z;
  t->bx_=(t->sumy_/t->nxy_)-t->ax_*(t->sumx_/t->nxy_);

     
  t->phi_=atan(t->ax_);
  if (t->phi_<0) t->phi_+=2*PI;

    double xp1=e->xp_[t->idx_[0]];
    double yp1=e->yp_[t->idx_[0]];
    if (xp1>0 && yp1>0 && t->phi_>PI) t->phi_-=PI;
    if (xp1<0 && yp1>0 && t->phi_>PI) t->phi_-=PI;
    if (xp1<0 && yp1<0 && t->phi_<PI) t->phi_+=PI;
    if (xp1>0 && yp1<0 && t->phi_<PI) t->phi_+=PI;

  t->theta_=atan(-1./t->ax_);
  t->pt_=5.7E-3*__fsqrt_ru((t->ax_*t->ax_+1)/t->bx_/t->bx_);
      //      printf("%f %f  \n",t->ax_,t->bx_);
  if (fabs(t->pt_)<1.8) {t->ok_=false;return;}
  s2z = t->sumz2_/t->nzr_-(t->sumz_/t->nzr_)*(t->sumz_/t->nzr_);
  szx = t->sumzr_/t->nzr_-(t->sumz_/t->nzr_)*(t->sumr_/t->nzr_);

  t->ar_= szx/s2z;
  t->br_=(t->sumr_/t->nzr_)-t->ar_*(t->sumz_/t->nzr_);
  t->z0_=-t->br_/t->ar_;
  if (fabs(t->z0_)>20) {t->ok_=false;return;}
  //t->eta_=-log(fabs(tan(atan( t->ar_)/2)));
  t->eta_=-log(fabs((1+__fsqrt_ru(1+t->ar_*t->ar_))/t->ar_) );
  if (t->ar_>0) t->eta_=-1*t->eta_;
  //if (t->nxy_>4) ngood++;
      
  t->chi2_=0;
  t->chi2r_=0;
  for (int ih=0;ih<t->nhit_;ih++)
    {
      int is=t->idx_[ih];
      double distx=(t->ax_*e->xp_[is]+t->bx_-e->yp_[is])/__fsqrt_ru(1+t->ax_*t->ax_)/8.4E-6;
      t->chi2_+=distx*distx;
      if (((e->lay_[is]>>16)&0x3)==0) continue;
      double deltar=(t->ar_*e->z_[is]+t->br_-e->r_[is])/__fsqrt_ru(1+t->ar_*t->ar_)/0.06;
      t->chi2r_+=deltar*deltar;

    }
    //__syncthreads();

    
}
void computeTklet(cudatklevent* e)
{
  //copyToHost(e);
  checkCudaErrors(cudaMemcpy(&e->host_->ntkl_,&e->device_->ntkl_,sizeof(unsigned int),cudaMemcpyDeviceToHost));
  dim3  grid2(e->host_->ntkl_,1, 1);
  //printf("%d %d %d\n",e->host_->ntkl_,e->host_->ntkl_/1024+1,min(e->host_->ntkl_,1024));
   dim3  grid1(e->host_->ntkl_/1024+1,1, 1);
   dim3  threads1(min(e->host_->ntkl_,1024), 1, 1);
  computeTkletKernel<<< grid1,threads1,0>>>(e->device_);
  //computeTkletKernel<<<grid2,1>>>(e->device_);
  cudaDeviceSynchronize();
    

}
