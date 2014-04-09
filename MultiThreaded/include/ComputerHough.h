#ifndef ComputerHough_H_
#define ComputerHough_H_
#include <vector>
#include <stdint.h>
#include "libhoughStruct.h"
#include "HoughStruct.h"
#include "libhoughCPU.h"
#include "TkComputer.h"



class ComputerHough : public virtual TkComputer
{
public:
  ComputerHough(HoughCut* cuts);
  virtual void DefaultCuts();
  virtual void Compute(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer);
  virtual void ComputeOneShot(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer);
  virtual void ComputeTracklet(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer,int32_t* flay);

  std::vector<mctrack_t> &getCandidates(){return theCandidateVector_;}
  static void Convert(double theta,double r,mctrack_t *m);
protected:
  std::vector<mctrack_t> theCandidateVector_;
  houghParam ph_;
  houghParam phcand_[96];
  houghParam phrcand_[64];
  HoughCut* theCuts_;
  uint32_t theNStub_;
  float* theX_;
  float* theY_;
  float* theZ_;
  uint32_t* theLayer_;
};
#endif
