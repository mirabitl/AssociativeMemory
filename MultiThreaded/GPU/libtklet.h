#ifndef _libtklet_h_
#define _libtklet_h_
#include "libhoughStruct.h"

typedef struct {

  unsigned int pattern_;
  unsigned int nhit_;
  unsigned short idx_[32];
  bool ok_,dummy_;
  double sumx_,sumx2_,sumxy_,sumy_,nxy_;
  double sumz_,sumz2_,sumzr_,sumr_,nzr_;
  double ax_,bx_,phi_,theta_,R_,pt_;
  double ar_,br_,eta_,z0_,chi2_,chi2r_;

} ctklet;

#define GPU_MAX_CTKLET 2048
typedef struct
{
  int flay_[32*128];
  float x_[GPU_MAX_STUB],y_[GPU_MAX_STUB];
  float z_[GPU_MAX_STUB],r_[GPU_MAX_STUB];
  float xp_[GPU_MAX_STUB],yp_[GPU_MAX_STUB];
  unsigned int lay_[GPU_MAX_STUB];
  unsigned int sector_,nstub_;
  unsigned int endcap_,barrel_,inter_;

  unsigned int ntkl_;
  ctklet cand_[GPU_MAX_CTKLET];

} ctklevent;

typedef struct 
{
  ctklevent* host_;
  ctklevent* device_;
} cudatklevent;

extern "C" void createTkletEvent(cudatklevent* e);
extern "C" void deleteTkletEvent(cudatklevent* e);
extern "C" void copyFromHost(cudatklevent* e);
extern "C" void copyToHost(cudatklevent* e);
extern "C" void fillDevice(cudatklevent* e);
extern "C" void combineLayer(cudatklevent* e);
extern "C" void addLayer(cudatklevent* e);
extern "C" void computeTklet(cudatklevent* e);
#endif
