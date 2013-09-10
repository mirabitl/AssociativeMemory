#ifndef _HOUGHLOCAL1_H

#define _HOUGHLOCAL1_H
#include <vector>
#include<stdint.h>
#include <DCHistogramHandler.h>
#define PI 3.141592653589793
#include "HoughStruct.h"
class RecoPoint;
class HoughVector 
{
 public:
  HoughVector(){this->clear();}
  void clear(){size_=0;memset(data_,0,32*sizeof(uint32_t));}
  void push_back(uint32_t i){if (size_<32) {data_[size_]=i;size_++;} }
  uint32_t size(){return size_;}
  inline uint32_t& operator[](uint32_t b){if (b<32) return data_[b];else {uint32_t err=0;return err;}}
 private:
  uint32_t size_,data_[32];
};
class HoughLocal1
{
public:
	
	HoughLocal1(double thmin,double thmax,double rmin,double rmax,uint32_t nbintheta=8,uint32_t nbinr=8);

	~HoughLocal1();
	void initialise(double thmin,double thmax,double rmin,double rmax,uint32_t nbintheta=8,uint32_t nbinr=8);
	void fill(double x,double y);
	void clear();
	double getTheta(int32_t i);
	double getR(int32_t i);
	uint16_t getValue(uint32_t i,uint32_t j);
	void findMaxima(std::vector< std::pair<uint32_t,uint32_t> >& maxval,uint32_t cut=3);
	void findMaxima(std::vector< std::pair<double,double> >& maxval,uint32_t cut=3);
	void findMaximumBins(std::vector< std::pair<double,double> >& maxval,uint32_t cut,std::vector<uint32_t> *ids=NULL);
	void findThresholdBins(std::vector< std::pair<double,double> >& maxval,uint32_t cut);

	void draw(DCHistogramHandler* h,std::vector< std::pair<uint32_t,uint32_t> > *maxval=NULL);
	void draw(DCHistogramHandler* h,std::vector< std::pair<double,double> > *maxval);
	uint32_t getVoteMax();
	static void PrintConvert(double theta,double r,mctrack_t *m=NULL);
	static void Convert(double theta,double r,mctrack_t *m);

	void addStub(stub_t s);
	// Inline getters
	inline double getThetaBin(){return  theThetaBin_;}
	inline double getRBin(){return theRBin_;}
	inline double getThetaMin(){return  theThetaMin_;}
	inline double getRMin(){return theRMin_;}
	inline double getThetaMax(){return  theThetaMax_;}
	inline double getRMax(){return theRMax_;}
	inline uint32_t getNbinTheta(){ return theNbinTheta_;}
	inline uint32_t getNbinR(){ return theNbinR_;}
	 
	inline HoughVector getHoughMap(uint16_t i,uint16_t j) { return theHoughMap_[i][j];}
	inline uint16_t getHoughImage(uint16_t i,uint16_t j) { return theHoughImage_[i][j];}
private:
	double theSin_[1024];
	double theCos_[1024];
	std::vector<double> theX_;
	std::vector<double> theY_;
	uint16_t theHoughImage_[1024][1024];
	HoughVector theHoughMap_[1024][1024];
	double theThetaMin_;
	double theThetaMax_;
	double theRMin_,theRMax_;
	double theThetaBin_,theRBin_;
	uint32_t theNbinTheta_;
	uint32_t theNbinR_;
	uint16_t theVoteMax_;
};



#endif
