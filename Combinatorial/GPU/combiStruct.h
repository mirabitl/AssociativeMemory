#ifndef _COMBISTRUCT_H
#define _COMBISTRUCT_H
#define GPU_MAX_STUB 128
#define GPU_MAX_LAYER 24
#define GPU_MAX_CAND 512
#ifndef PI
#define PI 3.141592653589793
#endif
typedef struct {
  float _x,_y,_z,_xp,_yp,_r2_r;
  unsigned int _id;
} stubPosition;
typedef struct {
  //Parameter
  int _nb;
  stubPosition stub[GPU_MAX_STUB];
} combiLayer;

static combiLayer *h_layer;
static combiLayer *d_layer;

typedef struct {
  double _x,_x2,_xy,_y,_nx;
  double _z,_z2,_zr,_r,_nz;
  double _a,_b,_ar,_br;
} combiTrack;

static combiTrack *h_cand;
static combiTrack *d_cand;




#endif
