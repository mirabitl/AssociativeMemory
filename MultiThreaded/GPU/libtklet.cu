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
  checkCudaErrors(cudaMemcpy(e->device_,e->host_,sizeof(ctklevent),
			     cudaMemcpyHostToDevice));
}

void copyToHost(cudatklevent* e)
{
  checkCudaErrors(cudaMemcpy(e->host_,e->device_,sizeof(ctklevent),
			     cudaMemcpyDeviceToHost));
}

__global__ void
fillConformalKernel(ctklevent* d)
{
  const unsigned int ib=blockIdx.x;

  float r2=d->x_[ib]*d->x_[ib]+d->y_[ib]*d->y_[ib];
  d->r_[ib]=__fsqrt_ru(r2);
  d->xp_[ib]=d->x_[ib]/r2;
  d->yp_[ib]=d->y_[ib]/r2;

}


void fillDevice(cudatklevent* e)
{
  copyFromHost(e);
  dim3  grid1(e->host_->nstub_, 1, 1);
  fillConformalKernel<<< grid1,1,0>>>(e->device_);
}

__global__ void
combineLayerKernel(ctklevent* d)
{
  const unsigned int is=blockIdx.x;
  const unsigned int js=threadIdx.x;
  const unsigned int l1=(d->lay_[is]&0xFFFF);
  const unsigned int l2=(d->lay_[js]&0xFFFF);
  bool barok= ((l1==5 && l2==6) || (l1==6 && l2==7) || (l1==5 && l2==7));
  if (d->barrel_ && !barok) return;
  if (d->inter_ && !barok) return;
  if (d->endcap_ && d->sector_<8)
    {
      if (!( (l1==5 && l2==18) || (l1==5 && l2==19) || (l1==18 && l2==19))) return;
    }
  if (d->endcap_ && d->sector_>=48)
    {
      if (!( (l1==5 && l2==11) || (l1==5 && l2==12) || (l1==11 && l2==12))) return;
    }

  
  // Ask z info
  if ( ((d->lay_[is]>>16)&0x3)==0) return;
  if ( ((d->lay_[js]>>16)&0x3)==0) return;
  //
  float a =(d->yp_[js]-d->yp_[is])/(d->xp_[js]-d->xp_[is]);
  //if ((isel%4==0) && (a<-0.25 || a>1.55)) continue;
  float b =d->yp_[js]-a*d->xp_[js];
  float pt=5.7E-3*__fsqrt_ru((a*a+1)/b/b);
  if (fabs(pt)<1.8) return;
  float ar=(d->r_[js]-d->r_[is])/(d->z_[js]-d->z_[is]);
  float br=d->r_[js]-ar*d->z_[js];
  float zi=-br/ar;
  if (fabs(zi)>20) return;
  unsigned int id=atomicInc(&d->ntkl_,GPU_MAX_CTKLET);
  ctklet* tkl=&d->cand_[id];
	  //printf("%x %x %d %d \n",d->lay_[i],d->lay_[j],(d->lay_[i]>>16)&0x3,(d->lay_[j]>>16)&0x3);
  tkl->nhit_=0; 
  tkl->pattern_=0; 
  tkl->ok_=true;
  tkl->pattern_ |=(1<<l1);
  tkl->pattern_ |=(1<<l2);
  tkl->idx_[tkl->nhit_++]=is;
  tkl->idx_[tkl->nhit_++]=js;

}
void combineLayer(cudatklevent* e)
{

  dim3  grid1(e->host_->nstub_,1, 1);
  combineLayerKernel<<< grid1,e->host_->nstub_,0>>>(e->device_);
}



__global__ void
addLayerKernel(ctklevent* e)
{
  const unsigned int is=blockIdx.x;
  const unsigned int it=threadIdx.x;
  ctklet* t=&(e->cand_[it]);
  if (!t->ok_) return;
  const unsigned int l1=(e->lay_[is]&0xFFFF);
  if ((t->pattern_&(1<<l1))!=0) return; // layer already included

  //double r2=e->r_[is]*e->r_[is];
  //r2=1.;
  double distx=(t->ax_*e->xp_[is]+t->bx_-e->yp_[is])/sqrt(1+t->ax_*t->ax_);
	  //	  if (fabs(distx)<1E-3)
	  //printf("is %d %f r2 %f %f %f %x pattern %x  \n",is,e->r_[is],distx,r2,t->ax_,e->lay_[is],t->pattern_);


	  //if (fabs(distx)>0.3) continue;
  
  double cut=3.5E-5;
  if (l1<=10) cut /=2;
  if (e->barrel_ && fabs(distx)>cut) return;
  if (e->inter_ && fabs(distx)>cut) return;
  if (e->endcap_ && fabs(distx)>cut) return;
  if (((e->lay_[is]>>16)&0x3)!=0)
    {
      double distr=(t->ar_*e->z_[is]+t->br_-e->r_[is])/sqrt(1+t->ar_*t->ar_);
      
      if (fabs(distr)>0.6) return;
	      //  if (e->barrel_ && fabs(distr)>0.45) continue;
	      // if (e->inter_ && fabs(distr)>0.2) continue;
	      //if (e->endcap_ && fabs(distr)>0.2) continue;
    }
  unsigned int id=atomicInc(&t->nhit_,32);
  t->idx_[id]=is;
  atomicOr(&t->pattern_,(1<<l1));
  t->pattern_|=(1<<l1);
}

void addLayer(cudatklevent* e)
{
  copyToHost(e);
  dim3  grid1(e->host_->nstub_,1, 1);
  dim3  grid2(e->host_->ntkl_,1, 1);
  addLayerKernel<<< grid1,grid2,0>>>(e->device_);
}

__global__ void
computeTkletKernel(ctklevent* e)
{
  const unsigned int it=threadIdx.x;
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
  float xp1=e->xp_[t->idx_[0]];
  float yp1=e->yp_[t->idx_[0]];
  if (xp1>0 && yp1>0 && t->phi_>PI) t->phi_-=PI;
  if (xp1<0 && yp1>0 && t->phi_>PI) t->phi_-=PI;
  if (xp1<0 && yp1<0 && t->phi_<PI) t->phi_+=PI;
  if (xp1>0 && yp1<0 && t->phi_<PI) t->phi_+=PI;
  t->theta_=atan(-1./t->ax_);
  t->pt_=5.7E-3*sqrt((t->ax_*t->ax_+1)/t->bx_/t->bx_);
      //      printf("%f %f  \n",t->ax_,t->bx_);
  if (fabs(t->pt_)<1.8) {t->ok_=false;return;}
  s2z = t->sumz2_/t->nzr_-(t->sumz_/t->nzr_)*(t->sumz_/t->nzr_);
  szx = t->sumzr_/t->nzr_-(t->sumz_/t->nzr_)*(t->sumr_/t->nzr_);

  t->ar_= szx/s2z;
  t->br_=(t->sumr_/t->nzr_)-t->ar_*(t->sumz_/t->nzr_);
  t->z0_=-t->br_/t->ar_;
  if (fabs(t->z0_)>20) {t->ok_=false;return;}
  //t->eta_=-log(fabs(tan(atan( t->ar_)/2)));
  t->eta_=-log(fabs((1+sqrt(1+t->ar_*t->ar_))/t->ar_) );
  if (t->ar_>0) t->eta_=-1*t->eta_;
  //if (t->nxy_>4) ngood++;
      
  t->chi2_=0;
  t->chi2r_=0;
  for (int ih=0;ih<t->nhit_;ih++)
    {
      int is=t->idx_[ih];
      double distx=(t->ax_*e->xp_[is]+t->bx_-e->yp_[is])/sqrt(1+t->ax_*t->ax_)/8.4E-6;
      t->chi2_+=distx*distx;
      if (((e->lay_[is]>>16)&0x3)==0) continue;
      double deltar=(t->ar_*e->z_[is]+t->br_-e->r_[is])/sqrt(1+t->ar_*t->ar_)/0.06;
      t->chi2r_+=deltar*deltar;

    }
    
    
}
void computeTklet(cudatklevent* e)
{
  copyToHost(e);

  dim3  grid2(e->host_->ntkl_,1, 1);
  computeTkletKernel<<< 1,grid2,0>>>(e->device_);
}
