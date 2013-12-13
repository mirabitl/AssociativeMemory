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



#include "combiStruct.h"

void allocateLayers()
{
  checkCudaErrors(cudaMalloc((void **) &d_layer,GPU_MAX_LAYER *sizeof(combiLayer)));
  checkCudaErrors(cudaMalloc((void **) &d_cand,GPU_MAX_CAND *sizeof(combiTrack)));
  checkCudaErrors(cudaHostAlloc((void **)&h_layer,GPU_MAX_LAYER *sizeof(combiLayer),cudaHostAllocWriteCombined));
  checkCudaErrors(cudaHostAlloc((void **)&h_cand,GPU_MAX_CAND *sizeof(combiTrack),cudaHostAllocWriteCombined));

}

void freeLayers()
{
  checkCudaErrors(cudaFreeHost(h_layer));
  checkCudaErrors(cudaFreeHost(h_cand));
  checkCudaErrors(cudaFree(d_layer));
  checkCudaErrors(cudaFree(d_cand));
}

void clearLayers()
{
  memset(h_layers,0,GPU_MAX_LAYER *sizeof(combiLayer));
  checkCudaErrors(cudaMemset(d_layer,0,GPU_MAX_LAYER *sizeof(combiLayer)));
}
void clearCandidates()
{
  memset(h_cand,0,GPU_MAX_CAND *sizeof(combiTrack));
  checkCudaErrors(cudaMemset(d_cand,0,GPU_MAX_CAND *sizeof(combiTrack)));
}

void copyLayer(uint32_t idl)
{
  checkCudaErrors(cudaMemcpy(&d_layers[idl], &h_layer[idl],
			     sizeof(int)+h_layer[idl]*sizeof(stubPosition),cudaMemcpyHostToDevice));
}

__global__ void
computeLayerKernel(combiLayer* L)
{
  const unsigned int ib=threadIdx.x;
  L->stub[ib]._r2= L->stub[ib]._x*L->stub[ib]._x+L->stub[ib]._y*L->stub[ib]._y;
  L->stub[ib]._r=sqrt(L->stub[ib]._r2);
  L->stub[ib]._xp=L->stub[ib]._x/L->stub[ib]._r2;
  L->stub[ib]._yp=L->stub[ib]._y/L->stub[ib]._r2;
  
}

void computeLayer(uint32_t idl)
{
  computeLayerKernel<<<1,h_layer[idl]._nb,0>>>(&d_layer[idl]);
}
