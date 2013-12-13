#include <iostream>
#include "TROOT.h"
#include "TSystem.h"

#include "L1TrackTrigger.h"
#include <boost/thread.hpp>  

#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TBrowser.h"
//TApplication* ta=NULL;
std::string filename;
void workerFunc(TApplication* ta)  
{  
  //boost::posix_time::seconds workTime(3);  
       
  std::cout << "Worker: running" << std::endl;  
  // Pretend to do something useful...  
  sleep((unsigned int) 1);
  ta->Run();
  //do not exit
  while (true)
    {
      //gSystem->Run();
					
      sleep((unsigned int) 30000000);
    }
  std::cout << "Worker: finished" << std::endl;  
} 

#undef FULLTREE 
void analysis(TApplication* theApp)
{
  //TBrowser b;
  //boost::thread workerThread(workerFunc,theApp);
	
#ifdef FULLTREE
  //getchar();
  TFile * datafile = new TFile("/home/mirabito/jets.root");
  TTree * awesomedata = (TTree *)datafile->Get("L1TrackTrigger");
  L1TrackTrigger l(awesomedata);
  l.Loop();
#else
  //getchar();
  L1TrackTrigger l;
  //std::string fname="/home/mirabito/activePatternsJetsBarrel64.root";
  //std::string fname="/home/mirabito/4on5/activePatternsBarrel64.root";
  //std::string fname="/home/mirabito/patterns6L/activePatterns140PU_4tops_6L_32ss_73cov_6on6.root";
  //std::string fname="/home/mirabito/patterns6LNoZ_16ss_30_6on6.root";
  //std::string fname="/home/mirabito/activePatternsBARREL_THRESH3_PU_5on6NUM.root";
  //std::string fname="/home/mirabito/activePatternsBARREL_THRESH2_PU_5on6_hough.root";
  std::string fname="/home/laurent/AM_Data/PU2_612_SLHC6_MUBANK_lowmidhig_sec40_ss32_cov40_4on6.root";
  std::cout<<filename<<std::endl;

#endif

}
int main(int argc, char* argv[])
{
  if (argc<3) exit(-1);
  std::cout<<"Directory "<<argv[1]<<std::endl;
  std::cout<<"File Name "<<argv[2]<<std::endl;
  std::string sdir(argv[1]);
  std::string sfile(argv[2]);
  filename=sdir+"/"+sfile+".root";
  std::cout<<filename<<std::endl;

  //ta=new 
  TApplication ta("THETEST",&argc,argv);
  //

  L1TrackTrigger l(0);
  l.do_ana(filename,0,sfile);


  //analysis(&ta);
  //workerThread.join();  
}
