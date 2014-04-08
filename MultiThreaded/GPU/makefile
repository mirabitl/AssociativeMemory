CUDA_EXAMPLE=${HOME}/cuda_examples/NVIDIA_CUDA-5.5_Samples/

libhough.o: libhough.cu libhough.h libhoughStruct.h
	/usr/local/cuda-5.5/bin/nvcc -ccbin g++ -I $(CUDA_EXAMPLE)/common/inc  -m64 -arch compute_30    -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=\"sm_35,compute_35\"  -o libhough.o -c libhough.cu

libtklet.o:libtklet.cu libtklet.h libhoughStruct.h
	/usr/local/cuda-5.5/bin/nvcc -ccbin g++ -I $(CUDA_EXAMPLE)/common/inc  -m64 -arch compute_30    -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=\"sm_35,compute_35\"  -o libtklet.o -c libtklet.cu
libhoughCPU.o: libhoughCPU.c libhoughCPU.h libhoughStruct.h
	g++ -c libhoughCPU.c
all:libhough.o libhoughCPU.o libtklet.o
