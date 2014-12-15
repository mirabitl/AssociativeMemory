#ifndef GenericAnalysis_H_
#define GenericAnalysis_H_
#include "DCHistogramHandler.h"
#include "HoughStruct.h"
#include <map>
#include <vector>
#include <list>
#undef USE_HV
#ifndef USE_HV
#define HOUGHLOCAL HoughLocal
#include "HoughLocal.h"
#else
#define HOUGHLOCAL HoughLocal1
#include "HoughLocal1.h"
#endif
#define WARN_PRINT_ENABLED 
#if DEBUG_PRINT_ENABLED
#define INFO_PRINT_ENABLED 1
#define DEBUG_PRINT fprintf
#else
#define DEBUG_PRINT(format, args...) ((void)0)
#endif
#if INFO_PRINT_ENABLED
#define WARN_PRINT_ENABLED 1
#define INFO_PRINT fprintf
#else
#define INFO_PRINT(format, args...) ((void)0)
#endif
#define WARN_PRINT_ENABLED 1
#if WARN_PRINT_ENABLED
#define WARN_PRINT printf
#else
#define WARN_PRINT(format, args...) ((void)0)
#endif
#ifdef USE_CUDA
#include "libhough.h"
#endif
#ifdef USE_CPU
#include "libhoughCPU.h"
#endif
typedef struct
{
  uint32_t goodmc;
  uint32_t missed;
  uint32_t fake;
  uint32_t rec;
} sectinfo;



class GenericAnalysis
{
public:
  enum FileType {GUILLAUME=1,SEBASTIEN=2};
  GenericAnalysis();
  void AddFile(std::string name,GenericAnalysis::FileType type);
  void FillMapGuillaumeNtuple(std::string name);
  void FillMapSebastienNtuple(std::string name);
  void CPULoopTest(std::string fname);
  void GBLoopTest(std::string fname);
  void ReadRawL1TrackTrigger(std::string fname,int evtshift=0);
  void ReadFullInfo(std::string fname,int mode=0);
  void MemoryLoopTest(std::string directory,int32_t sector=-1);
  void fill_histos();
  void associate();
  void event_hough(int isel);
  void basicHistos(int32_t i);
  void analyzePrecise();
  void alternativeAssociate();
  void PrintSectorMap();
  void cleanDuplicate(std::vector<mctrack_t> v);
#if defined(USE_CUDA) || defined(USE_CPU)
  void drawph(houghParam* p,DCHistogramHandler* r);
#endif
#if defined(USE_CUDA)
  void FillMapSectorNtuple(std::string fname);
  void FillMapOneShot(std::string fname);
  void FillMapEightSector(std::string fname);
#endif
protected:
  std::map<uint32_t,stub_t> theStubMap_;
  std::map<int32_t,mctrack_t> theMCMap_;
  std::vector<pattern_t> thePatternVector_;
  std::vector<mctrack_t> theHoughCandidateVector_;
  std::list<mctrack_t> theAssociatedTracks_;
  std::list<mctrack_t> theFakeTracks_;
  
  DCHistogramHandler* theRootHandler_;
  HOUGHLOCAL* theHoughLow_;
  HOUGHLOCAL* theHoughPrecise_;
  HOUGHLOCAL* theHoughR_;
  uint32_t ngoodmc_,nmiss_,nfake_;
  float thePtCut_,theNDelta_;
  uint32_t theNBinRho_;
  FILE* logFile_;
  bool barrel,endcap,inter;
  std::map<uint32_t,sectinfo> sectmap_;
  uint32_t theSector_;
};
#endif
