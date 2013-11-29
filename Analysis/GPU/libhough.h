#ifndef _LIBHOUGH_H
#define _LIBHOUGH_H
#define GPU_MAX_STUB 4096
#define GPU_MAX_THETA 256
#define GPU_MAX_RHO 256
#define GPU_MAX_RHO_WORD 16
#define GPU_MAX_STUB_BIN 16
#define GPU_MAX_CAND (GPU_MAX_THETA*GPU_MAX_RHO)
#define GPU_MAX_REG 100 
#define PI 3.141592653589793
#define RMIN -21.05
#define RMAX  21.05

typedef struct {
  //Parameter
  int nstub;
  int ntheta;
  int nrho;
  float thetamin;
  float thetamax;
  float thetabin;
  float rmin;
  float rmax;
  float rbin;
  // Device position
  float* d_x;
  float* d_y;
  float* d_z;
  float* d_r;
  unsigned int* d_layer;
  // Device image
  unsigned int* d_images;
  short* d_val;
  short* h_val;
  // Device hough
  unsigned int* d_hough;
  unsigned int* d_hough_map;
  unsigned int* d_hough_layer;

  unsigned int* d_temp;
  unsigned int* h_temp;
  // device points
  unsigned int* d_cand;
  // Host points
  unsigned int* h_cand;
  
  // Max value
  unsigned int max_val;

  // Rgression
  float* d_reg;
  float* h_reg;
  
} houghParam;


#define GET_R_VALUE(p,ir) (p.rmin+(ir+0.5)*p.rbin)
#define GET_THETA_VALUE(p,ith) (p.thetamin+(ith+0.5)*p.thetabin)

extern "C" void createHough(houghParam* p);
extern "C" void initialiseHough(houghParam* p,int nstub,int ntheta,int nrho,float thetamin,float thetamax,float rmin,float rmax);
extern "C" void deleteHough(houghParam* p);
extern "C" void fillPositionHough(houghParam* p,float* h_x,float* h_y,float* h_z);
extern "C" void fillConformalHough(houghParam* p,float* h_x,float* h_y,float* h_z);
extern "C" void fillLayerHough(houghParam* p,unsigned int* h_layer);

extern "C" void copyPositionHough(houghParam* pi,int icand,houghParam* po,unsigned int mode,bool regression);
extern "C" void doRegression(houghParam* p,unsigned int mode=0);
extern "C" void processHough(houghParam* p,unsigned int min_cut,unsigned int min_layer,unsigned int mode=0);
extern "C" void copyHoughImage(houghParam* p,unsigned int* h_hough);
extern "C" void copyHoughLayer(houghParam* p,unsigned int* h_hough);

extern "C" void clean(houghParam* p);
extern "C" void clearHough(houghParam* p);
extern "C" void dump(houghParam* p);
extern "C"  void initialiseTimer();
extern "C" void startTimer();
extern "C" float stopTimer();
extern "C" void deleteTimer();

#endif
