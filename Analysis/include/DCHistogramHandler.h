#ifndef _DCHISTOGRAMHANDLER_H
#define _DCHISTOGRAMHANDLER_H
#include <map>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <sys/types.h>
#include <regex.h>
#include <TMath.h>
#include "AbsTree.h"
#include <semaphore.h>
/** 
\class DCHistogramHandler
  \author  L.Mirabito 
  \date March 2010
  \version 1.0

   \brief ROOT online histogram handling 

   <h2> Description</h2>

   The <i>DCHistogramHandler</i> offers to book ROOT histograms in memory,
   keep  maps of histograms and write histograms to file in an online 
   environnement
*/
class DCHistogramHandler
{
 public:
  DCHistogramHandler();

  //! Book a TH1D histogram and register it in the map
  /**
    
     @param name string, it's the name of the histogram , its title and its entry in the map
     @param nbinx number of bin in x
     @param xmin low bin value
     @param xmax high bin value
     @return a TH1* ptr
     
   */
  TH1* BookTH1(std::string name,int nbinx,double xmin,double xmax); 
  //! Book a TH2D histogram and register it in the map
 /**
     
     @param name string, it's the name of the histogram , its title and its entry in the map
     @param nbinx number of bin in x
     @param xmin low bin value
     @param xmax high bin value
     @param nbiny number of bin in y
     @param ymin low y bin value
     @param ymax high y bin value
     @return a TH2* ptr
     
   */
  TH2* BookTH2(std::string name,int nbinx,double xmin,double xmax,int nbiny,double ymin,double ymax); 

  //!     Book a TProfile histogram and register it in the map (TH1 map)
 /**

     @param name string, it's the name of the histogram , its title and its entry in the map
     @param nbinx number of bin in x
     @param xmin low bin value
     @param xmax high bin value
     @param ymin min y value - if Ymax<Ymin no limit (default)
     @param ymax max y value
     @return a TProfile* ptr
     
   */
  TProfile* BookProfile(std::string name,int nbinx,double xmin,double xmax,double ymin=1,double ymax=-1); 

  //!     Book a TProfile2D histogram and register it in the map (TH2 map)
 /**

     @param name string, it's the name of the histogram , its title and its entry in the map
     @param nbinx number of bin in x
     @param xmin x low bin value
     @param xmax x high bin value
     @param nbiny number of bin in y
     @param ymin y low bin value 
     @param ymax y low bin value
     @param zmin min z value - if Zmax<Zmin no limit (default)
     @param zmax max z value
     @return a TProfile2D* ptr
     
   */

  TProfile2D* BookProfile2D(std::string name,int nbinx,double xmin,double xmax,int nbiny, double ymin,double ymax,double zmin=1,double zmax=-1); 

//! Book a TH3D histogram and register it in the map
 /**
     
     @param name string, it's the name of the histogram , its title and its entry in the map
     @param nbinx number of bin in x
     @param xmin low bin value
     @param xmax high bin value
     @param nbiny number of bin in y
     @param ymin low y bin value
     @param ymax high y bin value
     @param nbinz number of bin in z
     @param zmin low z bin value  - if Zmax<Zmin no limit
     @param zmax high z bin value
     @return a TH3* ptr
     
   */
  TH3* BookTH3(std::string name,int nbinx,double xmin,double xmax,int nbiny,double ymin,double ymax,int nbinz,double zmin,double zmax);

  //! Return a pointer to the TH1 or TProfile
  /** 
      @param name Name of the histogram
      @return a TH1* ptr
   */
  TH1* GetTH1(std::string name) 
{ if (mapH1.find(name)!=mapH1.end()) return mapH1[name]; else return 0;} 


  //! Return a pointer to the TH2
 /** 
      @param name Name of the histogram
      @return a Th2* ptr
   */
  TH2* GetTH2(std::string name) 
    { if (mapH2.find(name)!=mapH2.end()) return mapH2[name];else return 0; } 
  
  //! Return a pointer to the TH3
 /** 
      @param name Name of the histogram
      @return a Th3* ptr
   */
  TH3* GetTH3(std::string name) 
    { if (mapH3.find(name)!=mapH3.end()) return mapH3[name];else return 0; } 
  

  //! fill  a vector of names of histogram according to the regular experession
 /** 
      @param reg Linux regular expression
      @param type 1 or 2 according to histo type
      @param found reference to a vector of string
   */

  void ListHisto(std::string reg,int type,std::vector<std::string> &found); 

  //! Return  a TCanvas containing a drawing of the histograms whose name fit the regular expression
  /** 
      @param reg Linux regular expression
      @return a TCanvas* ptr
  */
  TCanvas* DrawRegexp(std::string reg,bool same=false); 

  //!Write histograms to a file
  /** 
      @param name File name
   */
  void writeHistograms(std::string name); 

  inline AbsTreeNode* getAbsTreeNode(){ return top;}

  //!Return XML list of Histo
  std::string getXMLHistoList();

  //!Return XML Buffer of a an Histo
  std::string getXMLHisto(std::string name);

  //!Lock and UnLock for multi threaded
  void Lock();
  void UnLock();
	//! Dump to XML
	void writeXML(std::string  path);
	//! referennce to the instance
  static DCHistogramHandler* instance() ;
 private:
  std::map<std::string,TH2*> mapH2;  /// Map of name , pointer to histogram TH1*
  std::map<std::string,TH1*> mapH1;  /// Map of name, pointer to histogram TH2*
  std::map<std::string,TH3*> mapH3;  /// Map of name, pointer to histogram TH3*

  AbsTreeNode* top;
  sem_t theMutex_;
  static DCHistogramHandler*  _me ;
};


#endif
