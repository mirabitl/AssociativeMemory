#ifndef _LIBHOUGH_H
#define _LIBHOUGH_H

#include "libhoughStruct.h"
extern "C" void createHough(houghParam* p,uint32_t max_stub=GPU_MAX_STUB,uint32_t max_theta=GPU_MAX_THETA,uint32_t max_rho=GPU_MAX_RHO);
extern "C" void initialiseHough(houghParam* p,int nstub,int ntheta,int nrho,float thetamin,float thetamax,float rmin,float rmax);
extern "C" void deleteHough(houghParam* p);
extern "C" void fillPositionHough(houghParam* p,float* h_x,float* h_y,float* h_z);
extern "C" void fillConformalHough(houghParam* p,float* h_x,float* h_y,float* h_z,int streamid=-1);
extern "C" void fillLayerHough(houghParam* p,unsigned int* h_layer,int streamid=-1);

extern "C" void copyPositionHough(houghParam* pi,int icand,houghParam* po,unsigned int mode,bool regression,int streamid=-1,bool endcap=false);
extern "C" void doRegression(houghParam* p,unsigned int mode=0);
extern "C" void processHough(houghParam* p,unsigned int min_cut,unsigned int min_layer,unsigned int mode,int streamid=-1,bool endcap=false);
extern "C" void copyHoughImage(houghParam* p,unsigned int* h_hough);
extern "C" void copyHoughLayer(houghParam* p,unsigned int* h_hough);

extern "C" void clean(houghParam* p);
extern "C" void getChi2(houghParam* p,bool endcap=false);
extern "C" void clearHough(houghParam* p);
extern "C" void dump(houghParam* p);
extern "C"  void initialiseTimer();
extern "C" void startTimer();
extern "C" float stopTimer();
extern "C" void deleteTimer();
extern "C" void* h_malloc(unsigned int memSize);
extern "C" void h_free(void* f);
extern "C" void createStreams(unsigned int nb);
extern "C" void deleteStreams(unsigned int nb);
extern "C" void  synchronize();

#endif
