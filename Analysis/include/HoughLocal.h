#ifndef _HOUGHLOCAL_H

#define _HOUGHLOCAL_H
#include <vector>
#include<stdint.h>
#include <DCHistogramHandler.h>
#define PI 3.141592653589793
#include "HoughStruct.h"
class RecoPoint;
class HoughLocal
{
public:
	
	HoughLocal(double thmin,double thmax,double rmin,double rmax,uint32_t nbintheta=8,uint32_t nbinr=8);

	~HoughLocal();
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
	
	inline std::vector<uint32_t> getHoughMap(uint16_t i,uint16_t j) { return theHoughMap_[i][j];}
	inline uint16_t getHoughImage(uint16_t i,uint16_t j) { return theHoughImage_[i][j];}
private:
	double theSin_[1024];
	double theCos_[1024];
	std::vector<double> theX_;
	std::vector<double> theY_;
	uint16_t theHoughImage_[1024][1024];
	std::vector<uint32_t> theHoughMap_[1024][1024];
	double theThetaMin_;
	double theThetaMax_;
	double theRMin_,theRMax_;
	double theThetaBin_,theRBin_;
	uint32_t theNbinTheta_;
	uint32_t theNbinR_;
	uint16_t theVoteMax_;
};



#endif
