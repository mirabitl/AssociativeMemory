#ifndef TkComputer_H_
#define TkComputer_H_

typedef struct
{
  float RhoMin;
  float RhoMax;
  uint32_t NRho;
  uint32_t NTheta;
  uint32_t NStubLow;
  uint32_t NLayerRow;
  uint32_t NStubLowCandidate;
  uint32_t NBins3GeV;
  uint32_t NBins5GeV; 
  uint32_t NBins15GeV;
  uint32_t NBins30GeV;
  uint32_t NBins100GeV;
  uint32_t NStubHigh;
  uint32_t NLayerHigh;
  uint32_t NStubHighCandidate;
  float NDelBarrel;
  float NDelInter;
  float NDelEndcap;
  
} HoughCut;



class TkComputer
{
public:
  virtual void DefaultCuts()=0;
  virtual void Compute(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer)=0;
  virtual void ComputeOneShot(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer)=0;
  virtual void ComputeTracklet(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer,int32_t* flay)=0;
  virtual void* getHough()=0;
};
#endif
