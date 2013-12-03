#include "GenericAnalysis.h"
#include <stdio.h>
#include <bitset>
#include <TROOT.h>
#include <TApplication.h>
#include <TChain.h>
#include <TFile.h>
#include <sstream>
#include <iostream>
using namespace std;
static TCanvas* CanvasGPU=NULL;
#ifdef USE_CUDA
#include "libhough.h"
#endif
#define USEMAIN
#ifdef USEMAIN
int main(int argc, char* argv[])
{
  TApplication ta("THETEST",&argc,argv);
  GenericAnalysis a;
  //a.AddFile("/home/mirabito/AM_Data/PU2_612_SLHC6_MUBANK_lowmidhig_sec16_ss32_cov40_5on6.root",GenericAnalysis::GUILLAUME);
  a.AddFile("/home/mirabito/AssociativeMemory/output_PU4TC_32_1000_COMPLETE.root",GenericAnalysis::SEBASTIEN);
  //a.AddFile("/home/mirabito/AssociativeMemory/output_PU_32_1000_ALL.root",GenericAnalysis::SEBASTIEN);
}
#endif
bool mctsort(mctrack_t t1,mctrack_t t2)
{
  if (t1.pt>t2.pt) return true;
  if (t1.pt == t2.pt) return (t1.phi>t2.phi);
  return false;
  
}
GenericAnalysis::GenericAnalysis() :thePtCut_(2.0),theNBinRho_(192),theNDelta_(1.5)
{
  theHoughLow_ = new HOUGHLOCAL(-PI/2,0.,-0.01,0.01,8,8);
  theHoughPrecise_ = new HOUGHLOCAL(-PI/2,0.,-0.01,0.01,8,8);
  theHoughR_ = new HOUGHLOCAL(-PI/2,0.,-0.01,0.01,8,8);
  theRootHandler_=DCHistogramHandler::instance();
  std::string sbase="/tmp/output";
  std::string logfilename="/tmp/output.log";
  logFile_ = fopen (logfilename.c_str() , "w");
}

void GenericAnalysis::AddFile(std::string name,GenericAnalysis::FileType type)
{
  if (type==GUILLAUME)
    FillMapGuillaumeNtuple(name);
  if (type==SEBASTIEN)
    FillMapSebastienNtuple(name);
}


void GenericAnalysis::FillMapSebastienNtuple(std::string fname)
{
  
  TChain *L1TT            = new TChain("FullInfo");  

  L1TT->Add(fname.c_str());

  //
  //1/ FullInfo TREE content:
  //
  // https://github.com/sviret/HL_LHC/blob/master/Extractors/RecoExtractor/test/SectorMaker/sector_test.h
  //


  int evt;             // Event number (for PU event, where there is more than 1 primary)
  int n_stub_total;    // The total number of stubs in the event
  int n_stub;          // The total number of stubs contained in matched patterns in the event

  std::vector<float>   *stub_x=new std::vector<float>;      // x coordinates of ALL the stubs
  std::vector<float>   *stub_y=new std::vector<float>;      // y coordinates of ALL the stubs
  std::vector<float>   *stub_z=new std::vector<float>;      // z coordinates of ALL the stubs
  std::vector<float>   *stub_x_2=new std::vector<float>;      // x coordinates of ALL the stubs
  std::vector<float>   *stub_y_2=new std::vector<float>;      // y coordinates of ALL the stubs
  std::vector<float>   *stub_z_2=new std::vector<float>;      // z coordinates of ALL the stubs
  std::vector<int>     *stub_layer=new std::vector<int>;  // layer number of ALL the stubs
  std::vector<int>     *stub_ladder=new std::vector<int>; // ladder number of ALL the stubs
  std::vector<int>     *stub_module=new std::vector<int>; // module number of ALL the stubs
  std::vector<int>     *stub_tp=new std::vector<int>;     // tp index of ALL the stubs (in part_*** vectors of this tree!!!!)
  std::vector<int>     *stub_inpatt=new std::vector<int>; // is the stub in a pattern (1) of not (0)?

  int n_part;                        // The total number of particles inducing at least one stub in the event

  std::vector<int>     *part_pdg=new std::vector<int>;    // PDG id of the particles
  std::vector<int>     *part_nsec=new std::vector<int>;   // In how many trigger towers this particle hit more than 4 different layers/disks?
  std::vector<int>     *part_nhits=new std::vector<int>;  // How many different layers/disks are hit by the particle?
  std::vector<int>     *part_npatt=new std::vector<int>;  // How many patterns contains more than 4 stubs of the particle (in 4 different layers/disks)?
  std::vector<float>   *part_pt=new std::vector<float>;     // pt of the particles
  std::vector<float>   *part_rho=new std::vector<float>;    // rho0 of the particles
  std::vector<float>   *part_z0=new std::vector<float>;     // z0 of the particles
  std::vector<float>   *part_eta=new std::vector<float>;    // eta of the particles 
  std::vector<float>   *part_phi=new std::vector<float>;    // phi of the particles

  int n_patt;                        // The total number of patterns matched in the event

  // Sector id of all the patterns
  std::vector<int>                  *patt_sec=new std::vector<int>;   

  // tp index of ALL the particles contained in the pattern (in part_*** vectors of this tree!!!!)
  std::vector< std::vector<int> >   *patt_parts=new std::vector< std::vector<int> >; 

  // index of ALL the stubs contained in the pattern (in stub_*** vectors of this tree!!!!) 
  std::vector< std::vector<int> >   *patt_stubs=new std::vector< std::vector<int> >; ; 

  L1TT->SetBranchAddress("evt",          &evt); 

  L1TT->SetBranchAddress("n_stub_total", &n_stub_total); 
  L1TT->SetBranchAddress("n_stub_inpat", &n_stub); 
  L1TT->SetBranchAddress("stub_x",       &stub_x); 
  L1TT->SetBranchAddress("stub_y",       &stub_y); 
  L1TT->SetBranchAddress("stub_z",       &stub_z); 
  L1TT->SetBranchAddress("stub_x_2",       &stub_x_2); 
  L1TT->SetBranchAddress("stub_y_2",       &stub_y_2); 
  L1TT->SetBranchAddress("stub_z_2",       &stub_z_2); 
  L1TT->SetBranchAddress("stub_layer",   &stub_layer); 
  L1TT->SetBranchAddress("stub_ladder",  &stub_ladder);
  L1TT->SetBranchAddress("stub_module",  &stub_module);
  L1TT->SetBranchAddress("stub_tp",      &stub_tp);
  L1TT->SetBranchAddress("stub_inpatt",  &stub_inpatt);

  L1TT->SetBranchAddress("n_part",       &n_part); 
  L1TT->SetBranchAddress("part_pdg",     &part_pdg); 
  L1TT->SetBranchAddress("part_nsec",    &part_nsec); 
  L1TT->SetBranchAddress("part_nhits",   &part_nhits); 
  L1TT->SetBranchAddress("part_npatt",   &part_npatt); 
  L1TT->SetBranchAddress("part_pt",      &part_pt); 
  L1TT->SetBranchAddress("part_rho",     &part_rho);
  L1TT->SetBranchAddress("part_z0",      &part_z0);
  L1TT->SetBranchAddress("part_eta",     &part_eta);  
  L1TT->SetBranchAddress("part_phi",     &part_phi); 

  L1TT->SetBranchAddress("n_patt",       &n_patt); 
  L1TT->SetBranchAddress("patt_sec",     &patt_sec); 
  L1TT->SetBranchAddress("patt_parts",   &patt_parts); 
  L1TT->SetBranchAddress("patt_stubs",   &patt_stubs); 

  int n_entries = L1TT->GetEntries();

  //if (evtnum>=n_entries) evtnum = n_entries-1;
  //if (evtnum<0) evtnum = 0;
  sectmap_.clear();
  for (uint32_t isect=1;isect<=56;isect++)
    {
      sectinfo s;
      memset(&s,0,sizeof(sectinfo));
      std::pair<uint32_t,sectinfo> p(isect,s);
      sectmap_.insert(p);
    }
  n_entries=1000;

#ifdef USE_CUDA  
  houghParam ph;
  createHough(&ph);
  //getchar();

  houghParam phi;
  createHough(&phi);
  // getchar();
  houghParam phreg;
  createHough(&phreg);

  houghParam phcand[64];
  for (int i=0;i<64;i++)
    createHough(&phcand[i]);

  createStreams(64);
  initialiseTimer();
#endif
  // allocate host memory useless but harmless for CPU
#ifdef USE_CUDA
  printf("On y va \n");
  float *h_x=(float *) h_malloc(1024*sizeof(float));
  float *h_y= (float *) h_malloc(1024*sizeof(float));
  float *h_z=(float *) h_malloc(1024*sizeof(float));
  unsigned int *h_layer=(unsigned int *) h_malloc(1024*sizeof(unsigned int));

  //  printf("Essai d'ecriture \n");
  // printf("Essai d'ecriture %f \n",h_x[129]);
  // getchar();
  // h_x[129]=12.;
  // printf("%f \n",h_x[129]);
  // getchar();
#else
  float *h_x=(float *) malloc(1024*sizeof(float));
  float *h_y= (float *) malloc(1024*sizeof(float));
  float *h_z=(float *) malloc(1024*sizeof(float));
  unsigned int *h_layer=(unsigned int *) malloc(1024*sizeof(unsigned int));

#endif
  float totalTime=0;

  for(int evtnum=1;evtnum<n_entries;evtnum++)
    {
      int gpu_nstub=0;
  
      L1TT->GetEntry(evtnum);

      cout <<endl;
      cout << "Analyzing event " << evt <<endl;
      cout << "where " << n_stub_total << " stub(s) were produced" <<endl;
      cout << n_part << " particle(s) have induced at least one stub in the tracker" <<endl;
      cout << "      " << n_stub << " stub(s) are contained in the " << n_patt << " pattern(s) matched in this event" <<endl;
      cout << endl;
      cout << "Now looping over the patterns... " <<endl;
      cout << endl;
      //getchar();
      // Loop over patterns

      int idx_s;
      int idx_p;
      for (int isel=16;isel<40;isel++)
	{
	  int gpu_nstub=0;
	  theHoughCandidateVector_.clear();
	  for (uint32_t isect=1;isect<57;isect++)
	    {
	      //if (isect!=9 && isect!=15 && isect!=43) continue;
	      //if (isect!=1  && isect!=5 && isect!=53) continue;
	      //if (isect!=16  && isect!=25 && isect!=38) continue;
	      if (isect!=isel) continue;
	  
	      endcap= (isect<16 || isect>=40);
	      barrel=!endcap;
	      if (endcap)
		{
		  theNBinRho_=160;
		  theNDelta_=2.5;
		}
	      else
		{
		  theNBinRho_=192;
		  theNDelta_=1.5;

		}
	      inter= (isect>=8 && isect<16) || (isect>=40 && isect<48);
	      if (inter) {barrel=true;endcap=false;}
	      theStubMap_.clear();
	      theMCMap_.clear();
	      theSector_=isect;
	      for (int k=0;k<n_patt;++k)
		{
		  if (patt_sec->at(k)!=isect ) continue;
		  /*
		    cout << "-------------------------------------------------"  <<endl;
		    cout << "Pattern " << k+1 << " properties:"  <<endl;
		    cout << "=> Sector id : " << patt_sec->at(k) <<endl;
		    cout << "=> Number of stubs : " << patt_stubs->at(k).size() <<endl;
		    cout << "=> Number of particles w/more than four stubs in the pattern : " << patt_parts->at(k).size() <<endl;
		  */
		  if (patt_parts->at(k).size()==0)
		    {
		      //cout << "!! FAKE PATTERN containing the following stubs: " <<endl;

		      for (int kk=0;kk<patt_stubs->at(k).size();++kk)
			{
			  idx_s = patt_stubs->at(k).at(kk);
			  idx_p = stub_tp->at(idx_s);
			  /*
			    cout << " Stub " << kk+1 << endl;  
			    cout << " X/Y/Z (in cm)       : " << stub_x->at(idx_s) 
			    << "/" << stub_y->at(idx_s) 
			    << "/" << stub_z->at(idx_s) << endl;
			  */

			  uint32_t hitIndex=idx_s;

			  //uint32_t sid=STUBID(stub_layer[hitIndex],stub_ladder[hitIndex],stub_z[hitIndex],stub_segment[hitIndex],stub_strip[hitIndex]);
			  std::map<uint32_t,stub_t>::iterator is=theStubMap_.find(idx_s);
			  if (is==theStubMap_.end())
			    {
	       
			      stub_t s;
			      s.id=idx_s;
			      s.x=stub_x->at(hitIndex);
			      s.y=stub_y->at(hitIndex);
			      s.z=stub_z->at(hitIndex);
			      s.r2=s.x*s.x+s.y*s.y;;
			      s.r=sqrt(s.r2);
			      s.xp=s.x/s.r2;
			      s.yp=s.y/s.r2;
			      s.tp=stub_tp->at(hitIndex);
			      s.layer =stub_layer->at(hitIndex);
			      //	DEBUG_PRINT(logFile_,"%d %d %f \n",theStubMap_.count(sid),hit_tp[hitIndex],hit_ptGEN[hitIndex]);
			      std::pair<uint32_t,stub_t> p(idx_s,s);
			      theStubMap_.insert(p);
			      h_x[gpu_nstub]=s.x;
			      h_y[gpu_nstub]=s.y;
			      h_z[gpu_nstub]=s.z;
			      h_layer[gpu_nstub]=s.layer;
			      gpu_nstub++;
#undef POINT2
#ifdef POINT2
			      h_x[gpu_nstub]=stub_x_2->at(hitIndex);
			      h_y[gpu_nstub]=stub_y_2->at(hitIndex);
			      h_z[gpu_nstub]=stub_z_2->at(hitIndex);
			      h_layer[gpu_nstub]=s.layer;
			      gpu_nstub++;
#endif
			    }
			  std::map<int32_t,mctrack_t>::iterator im=theMCMap_.find(stub_tp->at(hitIndex));
			  if (theMCMap_.count(stub_tp->at(hitIndex))==0 && stub_tp->at(hitIndex)>=0)
			    {
			      mctrack_t mct;
			      mct.id=stub_tp->at(hitIndex);
			      mct.phi=part_phi->at(mct.id);
			      if (mct.phi<0) mct.phi+=2*PI;
			      mct.pt=part_pt->at(mct.id);
			      mct.z0=part_z0->at(mct.id);
			      mct.rho0=part_rho->at(mct.id);
			      mct.eta=part_eta->at(mct.id);
			      mct.nhits=part_nhits->at(mct.id);
			      mct.nstubs=(1<<stub_layer->at(hitIndex));
			      mct.maxstubs=1;
			      mct.valid=false;
			      mct.matches=0;
			      // DEBUG_PRINT(logFile_,"insert done %d %d \n", mct.id,theMCMap_.count(mct.id));
			      std::pair<uint32_t,mctrack_t> mcp(mct.id,mct);
			      theMCMap_.insert(mcp);

			    }
			  else
			    if (im!=theMCMap_.end())
			      im->second.nstubs|=(1<<stub_layer->at(hitIndex));
	      
			  if (stub_tp->at(idx_s)>=0)  // The cluster is matched
			    {
			      /*
				cout << "    Matched with PART   : " << idx_p << endl;
				cout << "    PTgen     : " << part_pt->at(idx_p) 
				<< endl;
	  
				cout << "    PART origin (R/Z) : " << part_rho->at(idx_p) << " / " 
				<< part_z0->at(idx_p) << endl; 
				cout << "    PART pdg code       : " << part_pdg->at(idx_p) << endl; 
			      */
			    }
			  else
			    {
			      //cout << "Unmatched" << endl;
			    }
			}
		    }
		  else
		    {
		      for (int jj=0;jj<patt_parts->at(k).size();++jj)
			{
			  //cout << "Particle: " << jj+1 <<endl;
			}

		      /*
			cout << "This patterns contains the following stubs: " <<endl;
		      */
		      for (int kk=0;kk<patt_stubs->at(k).size();++kk)
			{
			  idx_s = patt_stubs->at(k).at(kk);
			  idx_p = stub_tp->at(idx_s);
			  /*
			    cout << " Stub " << kk+1 << endl;  
			    cout << " X/Y/Z (in cm)       : " << stub_x->at(idx_s) 
			    << "/" << stub_y->at(idx_s) 
			    << "/" << stub_z->at(idx_s) << endl;

			  */

			  uint32_t hitIndex=idx_s;

			  std::map<uint32_t,stub_t>::iterator is=theStubMap_.find(idx_s);
			  if (is==theStubMap_.end())
			    {
		      
			      stub_t s;
			      s.id=idx_s;

			      s.x=stub_x->at(hitIndex);
			      s.y=stub_y->at(hitIndex);
			      s.z=stub_z->at(hitIndex);
			      s.r2=s.x*s.x+s.y*s.y;;
			      s.r=sqrt(s.r2);
			      s.xp=s.x/s.r2;
			      s.yp=s.y/s.r2;
			      s.tp=stub_tp->at(hitIndex);
			      s.layer =stub_layer->at(hitIndex);
			      //	DEBUG_PRINT(logFile_,"%d %d %f \n",theStubMap_.count(sid),hit_tp[hitIndex],hit_ptGEN[hitIndex]);
			      std::pair<uint32_t,stub_t> p(idx_s,s);
			      theStubMap_.insert(p);
			      h_x[gpu_nstub]=s.x;
			      h_y[gpu_nstub]=s.y;
			      h_z[gpu_nstub]=s.z;
			      h_layer[gpu_nstub]=s.layer;
			      gpu_nstub++;
#ifdef POINT2
			      h_x[gpu_nstub]=stub_x_2->at(hitIndex);
			      h_y[gpu_nstub]=stub_y_2->at(hitIndex);
			      h_z[gpu_nstub]=stub_z_2->at(hitIndex);
			      h_layer[gpu_nstub]=s.layer;
			      gpu_nstub++;
#endif

			    }

			  std::map<int32_t,mctrack_t>::iterator im=theMCMap_.find(stub_tp->at(hitIndex));
			  if (theMCMap_.count(stub_tp->at(hitIndex))==0 && stub_tp->at(hitIndex)>=0)
			    {
			      mctrack_t mct;
			      mct.id=stub_tp->at(hitIndex);
			      mct.phi=part_phi->at(mct.id);
			      if (mct.phi<0) mct.phi+=2*PI;
			      mct.pt=part_pt->at(mct.id);
			      mct.z0=part_z0->at(mct.id);
			      mct.rho0=part_rho->at(mct.id);
			      mct.eta=part_eta->at(mct.id);
			      mct.nhits=part_nhits->at(mct.id);
			      mct.nstubs=(1<<stub_layer->at(hitIndex));
			      mct.maxstubs=1;
			      mct.valid=false;
			      mct.matches=0;
			      // DEBUG_PRINT(logFile_,"insert done %d %d \n", mct.id,theMCMap_.count(mct.id));
			      std::pair<uint32_t,mctrack_t> mcp(mct.id,mct);
			      theMCMap_.insert(mcp);

			    }
			  else
			    if (im!=theMCMap_.end())
			      im->second.nstubs|=(1<<stub_layer->at(hitIndex));




			  if (stub_tp->at(idx_s)>=0)  // The cluster is matched
			    {
			      /*
				cout << "    Matched with PART   : " << idx_p << endl;
				cout << "    PTgen     : " << part_pt->at(idx_p) 
				<< endl;
	  
				cout << "    PART origin (R/Z) : " << part_rho->at(idx_p) << " / " 
				<< part_z0->at(idx_p) << endl; 
				cout << "    PART pdg code       : " << part_pdg->at(idx_p) << endl; 
			      */
			    }
			  else
			    {
			      // cout << "Unmatched" << endl;
			    }
			}
		    }
		}
      
	      for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
		{

	    
		  uint8_t np=0;
		  for (uint8_t ib=5;ib<=24;ib++)
		    if ((im->second.nstubs>>ib)&1) np++;
		  //if (im->second.pt>5 ||im->second.pt<3 ) continue;	  
		  if (im->second.pt<thePtCut_-0.8) continue;	  
		  if (abs(im->second.rho0)>0.5) continue;	  
		  if (np>=5 &&im->second.pt>thePtCut_ && im->second.nhits>=5) 
		    {
		      INFO_PRINT("MC %d NSTUB %x %d PT %f   (phi) %f ->%d %d  %f\n",im->second.id,im->second.nstubs,np,im->second.pt,im->second.phi,np,im->second.nhits,im->second.z0);
		      if (np>im->second.maxstubs) 	im->second.maxstubs=np;
		      if (im->second.valid)					
			continue;
		      else
			{
			  im->second.valid=true;
			  //hext2->Fill(im->second.pt,tan(im->second.phi));
		
			  //hptgen->Fill(log(im->second.pt));
			  //hphigen->Fill(tan(im->second.phi));
			}
	    
		    }
		}

      
	      uint32_t ngood=0;
	      for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
		if (im->second.valid) ngood++;
	      DEBUG_PRINT(logFile_,"MC map size %d Good %d \n",(int) theMCMap_.size(),ngood);
	      INFO_PRINT("MC map size %d Good %d \n",(int) theMCMap_.size(),ngood);
	      //getchar();
	      // Make analysis
	      if (ngood==0) continue;

#ifdef USE_CUDA
	      startTimer();
	      if (gpu_nstub<1024) 
		{
		  float thmin=-PI/2,thmax=PI/2;
		  float rhmin=-0.0031,rhmax=0.0031;
		  INFO_PRINT("On appelle le GPU %d \n",gpu_nstub);
		  int ntheta=160;
		  int nrho=8;//12//192;
		  //initialiseHough(&ph,gpu_nstub,ntheta,nrho,-PI/2,PI/2,-0.06,0.06);
		  ntheta=48;//64;
		  if (isel%4==0) thmin=1.32;
		  if (isel%4==1) thmin=-1.04;
		  if (isel%4==2) thmin=-0.24;
		  if (isel%4==3) thmin=0.51;
		  thmax=thmin+1.25;
		  if (gpu_nstub>400)
		    {
		      ntheta*=2;
		      nrho*=2;
		    }
		  initialiseHough(&ph,gpu_nstub,ntheta,nrho,thmin,thmax,rhmin,rhmax);
		  //int nrho=12;
		  //initialiseHough(&ph,gpu_nstub,ntheta,nrho,-PI/2,PI/2,-0.004,0.004);

		  fillConformalHough(&ph,h_x,h_y,h_z);
		  fillLayerHough(&ph,h_layer);
		  //clearHough(&ph);
#ifndef POINT2
		  processHough(&ph,4,4,0,-1);//processHough(&ph,3,3,0);

#else
		  processHough(&ph,8,4,0);
#endif
		  INFO_PRINT("RMIN %f RMAX %f RBIN %f gives %d candidates Max val %d STubs %d\n",ph.rmin,ph.rmax,ph.rbin,ph.h_cand[0],ph.max_val,ph.nstub);
		  //if (ph.h_cand[0]>10) continue;
		  // if (ph.h_cand[0]>100)
		  //   {processHough(&ph,7);
		  //     printf("RMIN %f RMAX %f RBIN %f gives %d candidates \n",ph.rmin,ph.rmax,ph.rbin,ph.h_cand[0]);
		  //   }

		  //drawph(&ph, theRootHandler_);
		  TH1F* hdptreg=(TH1F*) theRootHandler_->GetTH1("dptreg");
		  TH1F* hdphireg=(TH1F*) theRootHandler_->GetTH1("dphireg");
		  TH1F* hnstubl=(TH1F*) theRootHandler_->GetTH1("nstubl");
		  TH1F* hnstubh=(TH1F*) theRootHandler_->GetTH1("nstubh");
		  TH1F* hncl=(TH1F*) theRootHandler_->GetTH1("ncl");
		  TH1F* hnch=(TH1F*) theRootHandler_->GetTH1("nch");
		  if (hdptreg==NULL)
		    {
		      hdptreg=(TH1F*)theRootHandler_->BookTH1("dptreg",200,-10.,10.);

		      hdphireg=(TH1F*)theRootHandler_->BookTH1("dphireg",200,-0.2,0.2);
		      hnstubl=(TH1F*)theRootHandler_->BookTH1("nstubl",500,0.,500.);
		      hnstubh=(TH1F*)theRootHandler_->BookTH1("nstubh",200,0.,200.);
		      hncl=(TH1F*)theRootHandler_->BookTH1("ncl",500,0.,500.);
		      hnch=(TH1F*)theRootHandler_->BookTH1("nch",200,0.,200.);
		    }
		  if (ph.h_cand[0]>0)
		    {
		      hnstubl->Fill(ph.nstub*1.);
		      hncl->Fill(ph.h_cand[0]*1.);
		    }
		  //doHough(nstub, h_x,h_y,ntheta,nrho,-PI/2.,PI/2,-21.,21.,h_cand);
		  /*
		    unsigned int candi[1024];
		    memcpy(candi,&ph.h_cand[1],ph.h_cand[0]*sizeof(unsigned int));
		    std::vector<unsigned int> vcand (candi, candi+ph.h_cand[0]);               // 32 71 12 45 26 80 53 33
		  */
		  // using default comparison (operator <):
		  //std::sort (vcand.begin(), vcand.end());


		      ///////////////////////////////////////////////////////
		  for (int ic=0;ic<TMath::Min(64,(int)ph.h_cand[0]);ic++)
		    {
		      phcand[ic].h_reg[20]=0;
		      int pattern=ph.h_cand[ic+1]; // vcand[ic]
		      int ith=pattern&0X3FF;
		      int ir=(pattern>>10)&0x3FF;
		      //ith=(vcand[ic])&0x3FF;
		      //ir=(vcand[ic]>>10)&0x3FF;
		      int ns=(pattern>>20)&0x3FF;
#ifndef POINT2
		      if (ns<5) continue;//if (ns<3) continue;
#else
		      if (ns<5) continue;
#endif
		      double PT=1./2./TMath::Abs(GET_R_VALUE(ph,ir))*0.3*3.8/100.;
		      if (PT<1.5) continue;
		      //printf("%f \n",TMath::Abs(GET_R_VALUE(ph,ir)));
		      uint32_t nbinf=64;
		      // <5,5-10,10-30,>30
		      if (PT<3) nbinf=56;
		      if (PT>=3 && PT<5) nbinf=128; // 128
		      if (PT>=5  && PT<15) nbinf=128;//192
		      if (PT>=15 && PT<=30) nbinf=128;//256
		      if (PT>=30 ) nbinf=128;//256

		      nbinf /=1;//1 //2 avant
		      //if (endcap) nbinf/=2;
		      uint32_t nbinr=nbinf;
#ifndef POINT2
		      if (ns>20 ) nbinf=2*nbinf;
		  
#else
		      if (ns>50 ) nbinf=2*nbinf;
#endif
		      float ndel=2.1;


		      float tmi=GET_THETA_VALUE(ph,ith)-ndel*ph.thetabin;

		      float tma=GET_THETA_VALUE(ph,ith)+ndel*ph.thetabin;
		      float rmi=GET_R_VALUE(ph,ir)-ndel*ph.rbin;
		      float rma=GET_R_VALUE(ph,ir)+ndel*ph.rbin;
	   
		      //printf(" From LowCandidat %f %d Look for bin  val= %x ns %d ith %d ir %d %f %f %f %f %d %d \n",PT,ph.max_val,pattern,ns,ith,ir,tmi,tma,rmi,rma,nbinf,nbinr);
		      //getchar();
		      nbinf/=1;
		      nbinr/=1;
		      //do {
			  initialiseHough(&phcand[ic],gpu_nstub,nbinf,nbinr,tmi,tma,rmi,rma);	    
			  //clearHough(&phi);
			  copyPositionHough(&ph,pattern,&phcand[ic],0,false,ic);
			  //printf(" From LowCandidat %f %d Look for bin  val= %x ns %d ith %d ir %d %f %f %f %f %d %d \n",PT,ph.max_val,pattern,phcand[ic].nstub,ith,ir,tmi,tma,rmi,rma,nbinf,nbinr);
		      
			  //dump(&phi);
		      //		getchar();
			/*
			unsigned int h_hough_l[ntheta*nrho];
			copyHoughLayer(&ph,h_hough_l);
			printf("Pattern %x \n",h_hough_l[ith*nrho+ir]);
			getchar();
			*/			
			  //printf("%d processed \n",ic);
		    }
		  
		  // getchar();
		  synchronize();
		  
		 
		  for (int ic=0;ic<TMath::Min(64,(int)ph.h_cand[0]);ic++)
		    {
		      if (phcand[ic].h_reg[20]>0)
			{
			  phcand[ic].nstub=int( phcand[ic].h_reg[20]);
			  processHough(&phcand[ic],5,5,0,ic);
			}
		    }
		  synchronize();

		  for (int ic=0;ic<TMath::Min(64,(int)ph.h_cand[0]);ic++)
		    {
		      if (phcand[ic].h_reg[20]>0)
			{
			  //printf("Candidats %d \n",phcand[ic].h_cand[0]);
			  //drawph(&phcand[ic],theRootHandler_);
			  // continue;
			  for (int ici=0;ici<phcand[ic].h_cand[0];ici++)
			    {
			      int patterni=phcand[ic].h_cand[ici+1]; 
			      int ithi=patterni&0X3FF;
			      int iri=(patterni>>10)&0x3FF;
			      
			      if (((patterni>>20)&0x3FF)<5) continue;
			      mctrack_t t;
			  
			      initialiseHough(&phreg,gpu_nstub,32,32,-PI/2,PI/2,-150.,150.);
			      copyPositionHough(&phcand[ic],patterni,&phreg,1,true);


			      if (phreg.h_reg[60+6]<1.7) continue;

			      if (phreg.h_reg[70+9]<1.5) continue;
			      t.z0=-phreg.h_reg[70+1]/phreg.h_reg[70+0];
			      t.eta=phreg.h_reg[70+8];
			      if (TMath::Abs(t.z0)>30.) continue;
			  

			      float theta=GET_THETA_VALUE(phcand[ic],ithi);
			      float r=GET_R_VALUE(phcand[ic],iri);

			      double a=-1./tan(theta);
			      double b=r/sin(theta);
			      
		
			  //
			      double R=1./2./TMath::Abs(r);
			      double xi=-a/2./b;
			      double yi=1./2./b;
			      double g_pt=0.3*3.8*R/100.;
			  //g_phi=atan(a);
			      double g_phi=theta-PI/2.;
			      if (g_phi<0) g_phi+=2*PI;
			      HOUGHLOCAL::Convert(theta,r,&t);
			      t.nhits=(patterni>>20)&0x3FF;
			  //t.z0=0;
			  theHoughCandidateVector_.push_back(t);

			    }
		      //doHough(nstub, h_x,h_y,64,64,tmi,tma,rmi,rma,h_cand1);











			}
		    }
		  goto endloop;
		      ///////////////////////////////////////////////////////


		  for (int ic=0;ic<TMath::Min(64,(int)ph.h_cand[0]);ic++)
		    {
		      int pattern=ph.h_cand[ic+1]; // vcand[ic]
		      int ith=pattern&0X3FF;
		      int ir=(pattern>>10)&0x3FF;
		      //ith=(vcand[ic])&0x3FF;
		      //ir=(vcand[ic]>>10)&0x3FF;
		      int ns=(pattern>>20)&0x3FF;
#ifndef POINT2
		      if (ns<5) continue;//if (ns<3) continue;
#else
		      if (ns<5) continue;
#endif
		      double PT=1./2./TMath::Abs(GET_R_VALUE(ph,ir))*0.3*3.8/100.;
		      if (PT<1.5) continue;
		      //printf("%f \n",TMath::Abs(GET_R_VALUE(ph,ir)));
		      uint32_t nbinf=64;
		      // <5,5-10,10-30,>30
		      if (PT<3) nbinf=56;
		      if (PT>=3 && PT<5) nbinf=128; // 128
		      if (PT>=5  && PT<15) nbinf=128;//192
		      if (PT>=15 && PT<=30) nbinf=128;//256
		      if (PT>=30 ) nbinf=128;//256

		      nbinf /=1;//1 //2 avant
		      //if (endcap) nbinf/=2;
		      uint32_t nbinr=nbinf;
#ifndef POINT2
		      if (ns>20 ) nbinf=2*nbinf;
		  
#else
		      if (ns>50 ) nbinf=2*nbinf;
#endif
		      float ndel=2.1;


		      float tmi=GET_THETA_VALUE(ph,ith)-ndel*ph.thetabin;

		      float tma=GET_THETA_VALUE(ph,ith)+ndel*ph.thetabin;
		      float rmi=GET_R_VALUE(ph,ir)-ndel*ph.rbin;
		      float rma=GET_R_VALUE(ph,ir)+ndel*ph.rbin;
	   
		      //printf(" From LowCandidat %f %d Look for bin  val= %x ns %d ith %d ir %d %f %f %f %f %d %d \n",PT,ph.max_val,pattern,ns,ith,ir,tmi,tma,rmi,rma,nbinf,nbinr);
		      //getchar();
		      nbinf/=1;
		      nbinr/=1;
		      //do {
			  initialiseHough(&phi,gpu_nstub,nbinf,nbinr,tmi,tma,rmi,rma);	    
			  //clearHough(&phi);
			  copyPositionHough(&ph,pattern,&phi,0,false);
			  // printf(" From LowCandidat %f %d Look for bin  val= %x ns %d ith %d ir %d %f %f %f %f %d %d \n",PT,ph.max_val,pattern,phi.nstub,ith,ir,tmi,tma,rmi,rma,nbinf,nbinr);
		      
			  //dump(&phi);
		      //		getchar();
			/*
			unsigned int h_hough_l[ntheta*nrho];
			copyHoughLayer(&ph,h_hough_l);
			printf("Pattern %x \n",h_hough_l[ith*nrho+ir]);
			getchar();
		      */
#ifndef POINT2
		      processHough(&phi,5,5,0);//processHough(&phi,4,4,0);
#else
		      processHough(&phi,8,5,0);
#endif
		      //drawph(&phi, theRootHandler_);
		      if (phi.h_cand[0]>30)
			{
			printf("%d %f GeV ns =%d HighPrec %d \n",ic,PT,phi.nstub,phi.h_cand[0]);
			//getchar();
			}
		      // if (phi.h_cand[0]>30)
		      // 	{
		      // 	  printf("%d ns =%d HighPrec %d \n",ic,phi.nstub,phi.h_cand[0]);
		      // 	  nbinf=int(nbinf*1.01);
		      // 	  nbinr=int(nbinr*1.01);
		      
		      // 	}
		      // }
		      // while (phi.h_cand[0]>50);


		      /*
			unsigned int icandi[1024];
			memcpy(icandi,&phi.h_cand[1],phi.h_cand[0]*sizeof(unsigned int));
			std::vector<unsigned int> vcandi (icandi, icandi+phi.h_cand[0]);               // 32 71 12 45 26 80 53 33

			// using default comparison (operator <):
			std::sort (vcandi.begin(), vcandi.end()); 
			//
			*/
		      if (phi.h_cand[0]>0)
			{
			  hnstubh->Fill(phi.nstub*1.);
			  hnch->Fill(phi.h_cand[0]*1.);
			}
		      for (int ici=0;ici<phi.h_cand[0];ici++)
			{
			  int patterni=phi.h_cand[ici+1]; 
			  int ithi=patterni&0X3FF;
			  int iri=(patterni>>10)&0x3FF;

			  if (((patterni>>20)&0x3FF)<5) continue;
			  //ithi=vcandi[ici]&0x3FF;
			  //iri=(vcandi[ici]>>10)&0x3FF;
#undef HOUGHZ
			  mctrack_t t;
			  
#define HOUGHZ_FIT
#ifdef HOUGHZ_FIT
			  initialiseHough(&phreg,gpu_nstub,32,32,-PI/2,PI/2,-150.,150.);
			  //clean(&phreg);
			  //clearHough(&phreg);
			  copyPositionHough(&phi,patterni,&phreg,1,true);
			  //dump(&phreg);
			  /*
			  doRegression(&phreg);
			  */

			  if (phreg.h_reg[60+6]<1.7) continue;
			  //regression Z

			  //doRegression(&phreg,1);

			  if (phreg.h_reg[70+9]<1.5) continue;
			  t.z0=-phreg.h_reg[70+1]/phreg.h_reg[70+0];
			  t.eta=phreg.h_reg[70+8];
			  if (TMath::Abs(t.z0)>30.) continue;
			  
#endif
			  float theta=GET_THETA_VALUE(phi,ithi);
			  float r=GET_R_VALUE(phi,iri);

			  double a=-1./tan(theta);
			  double b=r/sin(theta);
		
		
			  //
			  double R=1./2./TMath::Abs(r);
			  double xi=-a/2./b;
			  double yi=1./2./b;
			  double g_pt=0.3*3.8*R/100.;
			  //g_phi=atan(a);
			  double g_phi=theta-PI/2.;
			  if (g_phi<0) g_phi+=2*PI;
			  //#ifdef DRAW_DEBUG
			  //printf("ic %d ith %d ir %d %d From %d %d val=%x Pt %f \n",ic,ith,ir,ici,ithi,iri,vcandi[ici],g_pt);
		
			  if (g_pt>20000000.)
			    {
			      printf("%d From r=%f theta=%f a=%f b=%f  R= %f  => Pt=%f GeV/c  Phi0=%f \n",(patterni>>20)&0x3FF,r,theta,a,b,R,g_pt,g_phi);
			      for (int k=0;k<10;k++)
				printf("%f ",phreg.h_reg[k]);
			      printf("\n");
			      getchar();
			    }
			  //#endif
		

		
			  
	
			  HOUGHLOCAL::Convert(theta,r,&t);
			  //std::vector<uint32_t>::iterator itm=vn.begin();
			  //float x=theStubMap_[(*itm)].x;
			  //float y=theStubMap_[(*itm)].y;
			  //printf("%f %f %f \n",x,y,t.phi);
			  //if (x>0 && y>0 && t.phi>PI) t.phi-=PI;
			  //if (x<0 && y>0 && t.phi>PI) t.phi-=PI;
			  //if (x<0 && y<0 && t.phi<PI) t.phi+=PI;
			  //if (x>0 && y<0 && t.phi<PI) t.phi+=PI;


			  //t.phierr=theHoughPrecise_->getThetaBin();
			  //t.pterr=t.pt*theHoughPrecise_->getRBin()/abs(r);
			  //printf("%f +- %f  / %f +- %f \n",t.pt,t.pterr,t.phi,t.phierr);
			  /* if (x>0 && y<0 && t.phi<0) t.phi+=2*PI;
			     if (x>0 && y>0 && t.phi<0) t.phi+=2*PI;
			  */
#ifdef HOUGHZ
			  hdptreg->Fill(t.pt-phreg.h_reg[6]);
			  hdphireg->Fill(t.phi-phreg.h_reg[2]);
			  t.pt=phreg.h_reg[6];
			  t.phi=phreg.h_reg[2];
#endif
			  t.nhits=(patterni>>20)&0x3FF;
			  //t.z0=0;
#ifdef HOUGHZ
			  t.phierr=phi.thetabin;
			  r=phreg.h_reg[4];
			  t.r=r;
			  t.theta=theta;
			  t.pterr=t.pt*phi.rbin/abs(r);
#ifndef POINT2
			  processHough(&phreg,4,3,1);
#else
			  processHough(&phreg,8,3,1);
#endif
			  //printf("phreg %d %d \n",phreg.max_val,phreg.h_cand[0]);
			  //	getchar();
			  if (phreg.h_cand[0]==0) continue;
			  doRegression(&phreg,1);
#ifndef POINT2			  
			  if (phreg.h_reg[9]>=1.5) 
#else
			  if (phreg.h_reg[9]>=3.5) 
#endif
			    {
			      t.z0=-phreg.h_reg[1]/phreg.h_reg[0];
			      if (TMath::Abs(t.z0)>17.) continue;
			      //printf("phreg %d %f %f \n",phreg.nstub,phreg.h_reg[9],t.z0);
			    }
			  else
			    {
			      continue;
			      printf("phreg %d %f %f \n",phreg.nstub,phreg.h_reg[9],t.z0);
			      dump(&phreg);
			      getchar();
			    }
#endif
			  theHoughCandidateVector_.push_back(t);

			}
		      //doHough(nstub, h_x,h_y,64,64,tmi,tma,rmi,rma,h_cand1);

		      //getchar();
	   
		    }
		endloop:
		  INFO_PRINT("Fin du GPU %ld \n",	theHoughCandidateVector_.size() );

		}
	      //
	      //std::sort(theHoughCandidateVector_.begin(),theHoughCandidateVector_.end(),mctsort);
	      totalTime+=stopTimer();
#else
	      event_hough();
#endif

	      alternativeAssociate();
	      basicHistos(isel);
	      basicHistos(-1);
	      /*
	 
		INFO_PRINT("Fin du CPU\n");



	      */
	    }
	}
      if (evtnum%100 ==0)
	PrintSectorMap();
    }
  printf("TotalTime %f\n",totalTime);
  PrintSectorMap();
  std::string rfile="output_histos.root";
  theRootHandler_->writeHistograms(rfile);
#ifdef USE_CUDA
  deleteHough(&phi);
  deleteHough(&ph);
#endif
  //cudaDeviceReset();
}
void GenericAnalysis::FillMapGuillaumeNtuple(std::string fname)

{
  
   
   

  cout<<"On rentre"<<endl;
  cout<<"DCHist created"<<endl;
  
  const int MAX_NB_PATTERNS=1500;
  const int MAX_NB_HITS = 100;
  const int MAX_NB_LADDERS_PER_LAYER = 16;
  const int MAX_NB_LAYERS=6;

  // Il y a 2 TTree dans le fichier : le premier contient les secteurs avec un ID par secteur
  // le deuxième contient les evenements avec la liste des patterns par evenement ainsi que l'ID du secteur concerne
  TChain *PATT    = new TChain("Patterns"); // infos about patterns
  TChain *SEC    = new TChain("Sectors");   //infos about sectors

  PATT->Add(fname.c_str());
  SEC->Add(fname.c_str());

  int sector_id=0;
  int sector_layers=0;
  int nb_ladders_layer[MAX_NB_LAYERS];
  int sector_layer_list[MAX_NB_LAYERS];
  int sector_layer_0[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_1[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_2[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_3[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_4[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_5[MAX_NB_LADDERS_PER_LAYER];

  int* sector_layers_detail[MAX_NB_LAYERS];
  sector_layers_detail[0]=sector_layer_0;
  sector_layers_detail[1]=sector_layer_1;
  sector_layers_detail[2]=sector_layer_2;
  sector_layers_detail[3]=sector_layer_3;
  sector_layers_detail[4]=sector_layer_4;
  sector_layers_detail[5]=sector_layer_5;

  SEC->SetBranchAddress("sectorID",            &sector_id);//ID du secteur
  SEC->SetBranchAddress("nbLayers",            &sector_layers);// nombre de layers dans le secteur
  SEC->SetBranchAddress("layer",                  sector_layer_list);
  SEC->SetBranchAddress("nb_ladders_layer",       nb_ladders_layer); // nombre de ladders pour chaque layer
  SEC->SetBranchAddress("sectorLadders_layer_0",  sector_layer_0);//liste des ladders pour layer 0
  SEC->SetBranchAddress("sectorLadders_layer_1",  sector_layer_1);
  SEC->SetBranchAddress("sectorLadders_layer_2",  sector_layer_2);
  SEC->SetBranchAddress("sectorLadders_layer_3",  sector_layer_3);
  SEC->SetBranchAddress("sectorLadders_layer_4",  sector_layer_4);
  SEC->SetBranchAddress("sectorLadders_layer_5",  sector_layer_5);
	

  int nb_layers;
  int nb_patterns=0;
  int nb_tracks=0;
  int event_id;
  int superStrip_layer_0[MAX_NB_PATTERNS];
  int superStrip_layer_1[MAX_NB_PATTERNS];
  int superStrip_layer_2[MAX_NB_PATTERNS];
  int superStrip_layer_3[MAX_NB_PATTERNS];
  int superStrip_layer_4[MAX_NB_PATTERNS];
  int superStrip_layer_5[MAX_NB_PATTERNS];
  int pattern_sector_id[MAX_NB_PATTERNS];

  float track_pt[MAX_NB_PATTERNS];
  float track_phi[MAX_NB_PATTERNS];
  float track_d0[MAX_NB_PATTERNS];
  float track_eta[MAX_NB_PATTERNS];
  float track_z0[MAX_NB_PATTERNS];

  int nbHitPerPattern[MAX_NB_PATTERNS];
  short hit_layer[MAX_NB_PATTERNS*MAX_NB_HITS];
  short hit_ladder[MAX_NB_PATTERNS*MAX_NB_HITS];
  short hit_zPos[MAX_NB_PATTERNS*MAX_NB_HITS];
  short hit_segment[MAX_NB_PATTERNS*MAX_NB_HITS];
  short hit_strip[MAX_NB_PATTERNS*MAX_NB_HITS];
  int32_t hit_tp[MAX_NB_PATTERNS*MAX_NB_HITS];
  float hit_ptGEN[MAX_NB_PATTERNS*MAX_NB_HITS];
  float hit_etaGEN[MAX_NB_PATTERNS*MAX_NB_HITS];
  float hit_phiGEN[MAX_NB_PATTERNS*MAX_NB_HITS];
  float hit_x[MAX_NB_PATTERNS*MAX_NB_HITS];
  float hit_y[MAX_NB_PATTERNS*MAX_NB_HITS];
  float hit_z[MAX_NB_PATTERNS*MAX_NB_HITS];
  float hit_X0[MAX_NB_PATTERNS*MAX_NB_HITS];
  float hit_Y0[MAX_NB_PATTERNS*MAX_NB_HITS];
  float hit_Z0[MAX_NB_PATTERNS*MAX_NB_HITS];
	
  PATT->SetBranchAddress("nbLayers",            &nb_layers);//nombre de layers pour les patterns
  PATT->SetBranchAddress("nbPatterns",          &nb_patterns); // nombre de patterns dans l'evenement
  PATT->SetBranchAddress("nbTracks",            &nb_tracks);
  PATT->SetBranchAddress("eventID",             &event_id); // ID de l'evenement (le meme que dans le fichier de simulation)
  PATT->SetBranchAddress("sectorID",            pattern_sector_id);// ID du secteur du pattern (permet de retrouver le secteur dans le premier TTree)
  PATT->SetBranchAddress("superStrip0",         superStrip_layer_0);// tableau de superstrips pour le layer 0
  PATT->SetBranchAddress("superStrip1",         superStrip_layer_1);
  PATT->SetBranchAddress("superStrip2",         superStrip_layer_2);
  PATT->SetBranchAddress("superStrip3",         superStrip_layer_3);
  PATT->SetBranchAddress("superStrip4",         superStrip_layer_4);
  PATT->SetBranchAddress("superStrip5",         superStrip_layer_5);
  PATT->SetBranchAddress("nbStubs",             nbHitPerPattern); // nombre de stubs contenus dans chaque pattern
  PATT->SetBranchAddress("track_pt",            track_pt);//layer du stub
  PATT->SetBranchAddress("track_phi",           track_phi);//layer du stub
  PATT->SetBranchAddress("track_eta",           track_eta);//layer du stub
  PATT->SetBranchAddress("track_d0",            track_d0);//layer du stub
  PATT->SetBranchAddress("track_z0",            track_z0);//layer du stub
  PATT->SetBranchAddress("stub_layers",         hit_layer);//layer du stub
  PATT->SetBranchAddress("stub_ladders",        hit_ladder);//ladder du stub
  PATT->SetBranchAddress("stub_module",         hit_zPos);//position en Z du module du stub
  PATT->SetBranchAddress("stub_segment",        hit_segment);//segment du stub
  PATT->SetBranchAddress("stub_strip",          hit_strip);//numero de strip du stub
  PATT->SetBranchAddress("stub_tp",             hit_tp);//numero de la particule du stub
  PATT->SetBranchAddress("stub_ptGEN",          hit_ptGEN);//PT de la particule du stub
  PATT->SetBranchAddress("stub_etaGEN",         hit_etaGEN);//PT de la particule du stub
  PATT->SetBranchAddress("stub_phi0GEN",        hit_phiGEN);//PT de la particule du stub
  PATT->SetBranchAddress("stub_x",              hit_x);
  PATT->SetBranchAddress("stub_y",              hit_y);
  PATT->SetBranchAddress("stub_z",              hit_z);
  PATT->SetBranchAddress("stub_X0",             hit_X0);
  PATT->SetBranchAddress("stub_Y0",             hit_Y0);
  PATT->SetBranchAddress("stub_Z0",             hit_Z0);
	
  int* sstrips[6];
  sstrips[0]=superStrip_layer_0;
  sstrips[1]=superStrip_layer_1;
  sstrips[2]=superStrip_layer_2;
  sstrips[3]=superStrip_layer_3;
  sstrips[4]=superStrip_layer_4;
  sstrips[5]=superStrip_layer_5;


  /*
    Lecture du TTree des secteurs et chargement des infos dans une map
    sector_list[<ID du secteur>] retourne le secteur correspondant sous la forme d'un vector<vector<int>>
    c'est a dire un vector<int> par layer, le vector<int> contenant la liste des numeros de ladder
  */
  map<int,vector<vector<int> > > sector_list;
  int n_entries_SEC = SEC->GetEntries();
  cout<<"found "<<n_entries_SEC<<" entries for sectors"<<endl;
 



  for(int index=0;index<n_entries_SEC;index++){
    SEC->GetEntry(index);
    vector<vector<int> > s;
    for(int j=0;j<sector_layers;j++){
      vector<int> layer;
      for(int l=0;l<nb_ladders_layer[j];l++){
	layer.push_back(sector_layers_detail[j][l]);
      }
      s.push_back(layer);
    }
    sector_list[sector_id]=s;
  }
  delete SEC;
  ngoodmc_=0;nmiss_=0;nfake_=0;
  // Initi Histos
  TH2F*	hext2=(TH2F*) theRootHandler_->BookTH2("Ext2",200,0.,100.,200,-3.,3.);
  TH1F* hptgen=(TH1F*)theRootHandler_->BookTH1("PtGen",500,0.5,5.);
  TH1F* hphigen=(TH1F*)theRootHandler_->BookTH1("PhiGen",100,-3.,3.);
  TH1F* hncout=(TH1F*)theRootHandler_->BookTH1("ngen",20,0.,20.);
  TH1F*  hnfake=(TH1F*)theRootHandler_->BookTH1("nfake",1501,-0.9,1500.1);
  TH1F*  hntracks=(TH1F*)theRootHandler_->BookTH1("ntracks",51,-0.9,50.1);
  TH1F*  hnpatterns=(TH1F*)theRootHandler_->BookTH1("npatterns",1500,0.,1500.);
  TH1F*  hnmatches=(TH1F*)theRootHandler_->BookTH1("nmatches",500,-0.1,499.9);
  TH1F*  hnmisses=(TH1F*)theRootHandler_->BookTH1("nmisses",1500,0.,1500.);

  TH1F* hnstub=(TH1F*)theRootHandler_->BookTH1("nstub",200,0.,200.);

  TH1F* hpthough=(TH1F*)theRootHandler_->BookTH1("PtHough",500,0.5,5.);
  TH1F* hphihough=(TH1F*)theRootHandler_->BookTH1("PhiHough",100,-3.,3.);
  TH1F* hdpt=(TH1F*)theRootHandler_->BookTH1("deltaPt",200,-1.,1.);
  TH1F* hdphi=(TH1F*)theRootHandler_->BookTH1("deltaPhi",200,-0.05,0.05);
  TH2F*	hpt2=(TH2F*) theRootHandler_->BookTH2("Pt_Hough_Gen",200,0.,100.,200,0.,100.);
  TH2F*	hphi2=(TH2F*) theRootHandler_->BookTH2("Phi_Hough_Gen",200,-3.,3.,200,-3.,3.);
  TH2F*	hfake2=(TH2F*) theRootHandler_->BookTH2("Fake2",200,0.,100.,200,-3.,3.);
  TH2F*	hfake2e=(TH2F*) theRootHandler_->BookTH2("Fake2e",200,0.,1.,200,-3.,3.);
  
  TH2F*	hmiss2=(TH2F*) theRootHandler_->BookTH2("Miss2",200,0.,100.,200,-3.,3.);
  TH2F*	hfound2=(TH2F*) theRootHandler_->BookTH2("Found2",200,0.,100.,200,-3.,3.);
  TH2F*	hfound2e=(TH2F*) theRootHandler_->BookTH2("Found2e",200,-0.2,0.2,200,-0.05,0.05);
  TH1F*	herrors=(TH1F*)theRootHandler_->BookTH1("Errors",400,0.,8.E-3);	
  TProfile* hpt2p= theRootHandler_->BookProfile("Pt_Hough_Gen_Profile",200,0.,100.,0.,100.);
  TProfile* hphi2p=theRootHandler_->BookProfile("Phi_Hough_Gen_Profile",200,-3.,3.,-3.,3.);


  
  /*
    Lecture des patterns
  */
  int hitIndex = 0;
  int n_entries_MC = PATT->GetEntries();
  cout<<n_entries_MC<<" events found"<<endl<<endl;
  //HOUGHLOCAL* htl = new HOUGHLOCAL(0,PI,-0.01,0.01,64,64);
  HOUGHLOCAL* htl = new HOUGHLOCAL(-PI/2,0.,-0.01,0.01,96,96);
  //HoughCartesian* htl = new HoughCartesian(-0.0,0.03,0.0,0.03,128,32);
  //HoughRZ* htr = new HoughRZ(PI/2,3*PI/2,-20,20,512,512);
  //float xpos[1000],ypos[1000];
  // Boucle sur les evenements
  std::map<uint32_t,uint32_t> counts;

  uint32_t nevmax=0;
  if (nevmax!=0 && nevmax<n_entries_MC)
    n_entries_MC=nevmax;
  uint32_t nfaketot=0,nfaketot1=0;
  time_t t0=time(0);
  DEBUG_PRINT(logFile_,"The tree has %d entries \n",n_entries_MC);
  for (int j=0;j<n_entries_MC;++j)
    {
      PATT->GetEntry(j); // Load entries
      theHoughCandidateVector_.clear();

      if(nb_patterns>0){//On a au moins un pattern actif
	//if (nb_tracks!=1) continue;
	//if (track_pt[0]>25) continue;
	time_t t1=time(0);
	cout<<nb_patterns<<"New EVENT pattern(s) in evt "<<event_id<<" (index "<<j<<")"<<endl;
	if (j%5==0) std::cout<<t1-t0<<"s Rate "<<j*1./(t1-t0)<<endl;
	hnpatterns->Fill(nb_patterns*1.);
	hitIndex=0;
			
	//boucle sur les patterns actifs de l'evenement
	//memset(xpos,0,1000*sizeof(float));
	//memset(ypos,0,1000*sizeof(float));
	uint32_t hitf=0;
	uint32_t hitl=0;
	uint32_t nfake=0;
	// FInd all MC tracks with >=3 stubs in a given patterns
	theMCMap_.clear();
	theStubMap_.clear();
	thePatternVector_.clear();
	TH1F* hnfake=(TH1F*)theRootHandler_->GetTH1("nfake");
	TH2F* hmiss2=(TH2F*) theRootHandler_->GetTH2("Miss2");
	TH2F* hfound2=(TH2F*) theRootHandler_->GetTH2("Found2");
 
	for(int k=0;k<nb_patterns;k++){
	  //if (nbHitPerPattern[k]<5) continue;

	  // Clear stubs count

	  pattern_t pats;
	  for(int m=0;m<nbHitPerPattern[k];m++)
	    {
	      //cout<<"Layer "<<(int)hit_layer[hitIndex]<<" ladder "<<(int)hit_ladder[hitIndex]<<" Mod "<<(int)hit_zPos[hitIndex]
	      // <<" Segment "<<(int)hit_segment[hitIndex]<<" strip "<<hit_strip[hitIndex]<<" (tp : "<<hit_tp[hitIndex]<<" PT : "<<hit_ptGEN[hitIndex]<<" GeV ETA : "
	      // <<hit_etaGEN[hitIndex]<<" PHI : "<<hit_phiGEN[hitIndex]<<" X:"<<hit_x[hitIndex]<<" Y:"<<hit_y[hitIndex]<<" Z:"<<hit_z[hitIndex]<<")"<<endl;
	      uint32_t sid=STUBID(hit_layer[hitIndex],hit_ladder[hitIndex],hit_zPos[hitIndex],hit_segment[hitIndex],hit_strip[hitIndex]);
	      std::map<uint32_t,stub_t>::iterator is=theStubMap_.find(sid);
	      if (is==theStubMap_.end())
		{
	       
		  stub_t s;
		  s.id=sid;
		  s.x=hit_x[hitIndex];
		  s.y=hit_y[hitIndex];
		  s.z=hit_z[hitIndex];
		  s.r2=s.x*s.x+s.y*s.y;;
		  s.xp=s.x/s.r2;
		  s.yp=s.y/s.r2;
		  s.tp=hit_tp[hitIndex];
		  s.layer=hit_layer[hitIndex];
		  //	DEBUG_PRINT(logFile_,"%d %d %f \n",theStubMap_.count(sid),hit_tp[hitIndex],hit_ptGEN[hitIndex]);
		  std::pair<uint32_t,stub_t> p(sid,s);
		  theStubMap_.insert(p);
		}
	      std::map<int32_t,mctrack_t>::iterator im=theMCMap_.find(hit_tp[hitIndex]);
	      if (theMCMap_.count(hit_tp[hitIndex])==0 && hit_tp[hitIndex]>=0)
		{
		  mctrack_t mct;
		  mct.id=hit_tp[hitIndex];
		  mct.phi=hit_phiGEN[hitIndex];
		  //if (mct.phi<0) mct.phi+=2*PI;
		  mct.pt=hit_ptGEN[hitIndex];
		  mct.z0=hit_Z0[hitIndex];
		  mct.eta=hit_etaGEN[hitIndex];
		  mct.nstubs=(1<<hit_layer[hitIndex]);
		  mct.maxstubs=1;
		  mct.valid=false;
		  mct.matches=0;
		  // DEBUG_PRINT(logFile_,"insert done %d %d \n", mct.id,theMCMap_.count(mct.id));
		  std::pair<uint32_t,mctrack_t> mcp(mct.id,mct);
		  theMCMap_.insert(mcp);

		}
	      else
		if (im!=theMCMap_.end())
		  im->second.nstubs|=(1<<hit_layer[hitIndex]);
	      pats.stubs_id.push_back(sid);
	      hitIndex++;

	    }
	
	  //DEBUG_PRINT(logFile_,"Pattern %d Nb stubs %d Total stubs %d Total MC %d \n",k,pats.stubs_id.size(),theStubMap_.size(),theMCMap_.size());
	  thePatternVector_.push_back(pats);
	  //getchar();
	}

	//for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
	//	im->second.nstubs=0;
      
	for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
	  {


	    uint8_t np=0;
	    for (uint8_t ib=5;ib<=24;ib++)
	      if ((im->second.nstubs>>ib)&1) np++;
	
	    if (im->second.pt<thePtCut_-0.8) continue;	  
	    if (np>=4 &&im->second.pt>thePtCut_ ) 
	      {
		DEBUG_PRINT(logFile_,"MC %d NSTUB %x %d PT %f   tg(phi) %f\n",im->second.id,im->second.nstubs,np,im->second.pt,tan(im->second.phi));
		if (np>im->second.maxstubs) 	im->second.maxstubs=np;
		if (im->second.valid)					
		  continue;
		else
		  {
		    im->second.valid=true;
		    hext2->Fill(im->second.pt,tan(im->second.phi));
		
		    hptgen->Fill(log(im->second.pt));
		    hphigen->Fill(tan(im->second.phi));
		  }
	    
	      }
	  }

      
	uint32_t ngood=0;
	for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
	  if (im->second.valid) ngood++;
	DEBUG_PRINT(logFile_,"MC map size %d Good %d \n",(int) theMCMap_.size(),ngood);
	event_hough();
	associate();
	//fill_histos(&rootHandler_);
	continue;
      }
		
    }
  delete PATT;

  INFO_PRINT(logFile_,"Good MC= %d  Missed = %d Fake = %d \n",ngoodmc_,nmiss_,nfake_);

  std::string rfile="output_histos.root";
  theRootHandler_->writeHistograms(rfile);
}

void GenericAnalysis::basicHistos(int32_t isel)
{
  if (isel>=0)
    {
      std::stringstream s;
      s<<"/sector"<<isel<<"/";
  
      TH1F* hptgood=(TH1F*) theRootHandler_->GetTH1(s.str()+"Ptgood");
      TH1F* hptgoodfound=(TH1F*) theRootHandler_->GetTH1(s.str()+"Ptgoodfound");
      TH1F* hphigood=(TH1F*) theRootHandler_->GetTH1(s.str()+"Phigood");
      TH1F* hphigoodfound=(TH1F*) theRootHandler_->GetTH1(s.str()+"Phigoodfound");
      TH1F* hdpt=(TH1F*) theRootHandler_->GetTH1(s.str()+"dpt");
      TH1F* hdptrel=(TH1F*) theRootHandler_->GetTH1(s.str()+"dptrel");
      TH1F* hdphi=(TH1F*) theRootHandler_->GetTH1(s.str()+"dphi");
      TH1F* hdz=(TH1F*) theRootHandler_->GetTH1(s.str()+"dz");
      TH1F* hzt=(TH1F*) theRootHandler_->GetTH1(s.str()+"zt");
      TH1F* hr=(TH1F*) theRootHandler_->GetTH1(s.str()+"r");
      TH1F* htheta=(TH1F*) theRootHandler_->GetTH1(s.str()+"theta");
      TH1F* hzmc=(TH1F*) theRootHandler_->GetTH1(s.str()+"zmc");
      TH2F* hdptnear_pth = (TH2F*) theRootHandler_->GetTH2(s.str()+"Dptnear_pth");
      TH2F* hptgen_pth = (TH2F*) theRootHandler_->GetTH2(s.str()+"ptgen_pth");
      TH2F* hdptonear_pth = (TH2F*) theRootHandler_->GetTH2(s.str()+"Dptonear_pth");
      TH2F* hdphinear_phih = (TH2F*) theRootHandler_->GetTH2(s.str()+"Dphinear_phi_h");
      if (hptgood==NULL)
	{
	  hr=(TH1F*)theRootHandler_->BookTH1(s.str()+"r",200,-0.05,0.05);
	  htheta=(TH1F*)theRootHandler_->BookTH1(s.str()+"theta",200,-PI/2,PI/2);


	  hptgood=(TH1F*)theRootHandler_->BookTH1(s.str()+"Ptgood",150,0.,50.);
	  hptgoodfound=(TH1F*)theRootHandler_->BookTH1(s.str()+"Ptgoodfound",150,0.,50.);
	  hphigood=(TH1F*)theRootHandler_->BookTH1(s.str()+"Phigood",100,0.,2*PI);
	  hphigoodfound=(TH1F*)theRootHandler_->BookTH1(s.str()+"Phigoodfound",100,0.,2*PI);
	  hdpt=(TH1F*)theRootHandler_->BookTH1(s.str()+"dpt",200,-10.,10.);
	  hdptrel=(TH1F*)theRootHandler_->BookTH1(s.str()+"dptrel",200,-0.5,0.5);
	  hdphi=(TH1F*)theRootHandler_->BookTH1(s.str()+"dphi",200,-0.2,0.2);
	  hdz=(TH1F*)theRootHandler_->BookTH1(s.str()+"dz",200,-2.,2.);
	  hzmc=(TH1F*)theRootHandler_->BookTH1(s.str()+"zmc",200,-100.,100.);
	  hzt=(TH1F*)theRootHandler_->BookTH1(s.str()+"zt",200,-100.,100.);
	  hdptnear_pth= (TH2F*) theRootHandler_->BookTH2(s.str()+"Dptnear_pth",50,0.,50.,100,-10.,10.);
	  hptgen_pth= (TH2F*) theRootHandler_->BookTH2(s.str()+"ptgen_pth",50,0.,50.,50,0.,50.);
      
	}
  for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
    if (im->second.valid)
      {
	hzmc->Fill(im->second.z0);
	hptgood->Fill(im->second.pt);
	if (im->second.matches) hptgoodfound->Fill(im->second.pt);
	hphigood->Fill(im->second.phi);
	if (im->second.matches) hphigoodfound->Fill(im->second.phi);
      }
  for (std::vector<mctrack_t>::iterator it=theAssociatedTracks_.begin();it!=theAssociatedTracks_.end();it++)
    {
      hr->Fill(it->r);
      hzt->Fill(it->z0);
      htheta->Fill(it->theta);
      hptgen_pth->Fill(it->pt,theMCMap_[it->id_ass].pt);
      hdptnear_pth->Fill(it->pt,it->pt-theMCMap_[it->id_ass].pt);
      hdpt->Fill(it->pt-theMCMap_[it->id_ass].pt);
      hdptrel->Fill((it->pt-theMCMap_[it->id_ass].pt)/it->pt);
      hdphi->Fill(it->phi-theMCMap_[it->id_ass].phi);
      hdz->Fill(it->z0-theMCMap_[it->id_ass].z0);
    }

    }
  else
    {
      TH1F* hptgood=(TH1F*) theRootHandler_->GetTH1("Ptgood");
      TH1F* hptgoodfound=(TH1F*) theRootHandler_->GetTH1("Ptgoodfound");
      TH1F* hphigood=(TH1F*) theRootHandler_->GetTH1("Phigood");
      TH1F* hphigoodfound=(TH1F*) theRootHandler_->GetTH1("Phigoodfound");
      TH1F* hdpt=(TH1F*) theRootHandler_->GetTH1("dpt");
      TH1F* hdptrel=(TH1F*) theRootHandler_->GetTH1("dptrel");
      TH1F* hdphi=(TH1F*) theRootHandler_->GetTH1("dphi");
      TH1F* hdz=(TH1F*) theRootHandler_->GetTH1("dz");
      TH1F* hdeta=(TH1F*) theRootHandler_->GetTH1("deta");
      TH1F* hr=(TH1F*) theRootHandler_->GetTH1("r");
      TH1F* htheta=(TH1F*) theRootHandler_->GetTH1("theta");
      TH1F* hzmc=(TH1F*) theRootHandler_->GetTH1("zmc");
      TH1F* hzt=(TH1F*) theRootHandler_->GetTH1("zt");
      TH2F* hdptnear_pth = (TH2F*) theRootHandler_->GetTH2("Dptnear_pth");
      TH2F* hptgen_pth = (TH2F*) theRootHandler_->GetTH2("ptgen_pth");
      TH2F* hdptonear_pth = (TH2F*) theRootHandler_->GetTH2("Dptonear_pth");
      TH2F* hdphinear_phih = (TH2F*) theRootHandler_->GetTH2("Dphinear_phi_h");
      if (hptgood==NULL)
	{
	  hr=(TH1F*)theRootHandler_->BookTH1("r",200,-0.05,0.05);
	  htheta=(TH1F*)theRootHandler_->BookTH1("theta",200,-PI/2,PI/2);


	  hptgood=(TH1F*)theRootHandler_->BookTH1("Ptgood",150,0.,50.);
	  hptgoodfound=(TH1F*)theRootHandler_->BookTH1("Ptgoodfound",150,0.,50.);
	  hphigood=(TH1F*)theRootHandler_->BookTH1("Phigood",100,0.,2*PI);
	  hphigoodfound=(TH1F*)theRootHandler_->BookTH1("Phigoodfound",100,0.,2*PI);
	  hdpt=(TH1F*)theRootHandler_->BookTH1("dpt",200,-10.,10.);
	  hdptrel=(TH1F*)theRootHandler_->BookTH1("dptrel",200,-0.5,0.5);
	  hdphi=(TH1F*)theRootHandler_->BookTH1("dphi",200,-0.2,0.2);
	  hdz=(TH1F*)theRootHandler_->BookTH1("dz",200,-2.,2.);
	  hdeta=(TH1F*)theRootHandler_->BookTH1("deta",200,-0.75,0.75);
	  hzmc=(TH1F*)theRootHandler_->BookTH1("zmc",200,-100.,100.);
	  hzt=(TH1F*)theRootHandler_->BookTH1("zt",200,-100.,100.);
	  hdptnear_pth= (TH2F*) theRootHandler_->BookTH2("Dptnear_pth",50,0.,50.,100,-10.,10.);
	  hptgen_pth= (TH2F*) theRootHandler_->BookTH2("ptgen_pth",50,0.,50.,50,0.,50.);
      
	}
  for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
    if (im->second.valid)
      {
	hzmc->Fill(im->second.z0);
	hptgood->Fill(im->second.pt);
	if (im->second.matches) hptgoodfound->Fill(im->second.pt);
	hphigood->Fill(im->second.phi);
	if (im->second.matches) hphigoodfound->Fill(im->second.phi);
      }
  for (std::vector<mctrack_t>::iterator it=theAssociatedTracks_.begin();it!=theAssociatedTracks_.end();it++)
    {
      hr->Fill(it->r);
      hzt->Fill(it->z0);
      htheta->Fill(it->theta);
      hptgen_pth->Fill(it->pt,theMCMap_[it->id_ass].pt);
      hdptnear_pth->Fill(it->pt,it->pt-theMCMap_[it->id_ass].pt);
      hdpt->Fill(it->pt-theMCMap_[it->id_ass].pt);
      hdptrel->Fill((it->pt-theMCMap_[it->id_ass].pt)/it->pt);
      hdphi->Fill(it->phi-theMCMap_[it->id_ass].phi);
      hdz->Fill(it->z0-theMCMap_[it->id_ass].z0);
      hdeta->Fill(it->eta-theMCMap_[it->id_ass].eta);
      //printf("%f %f \n",it->eta,theMCMap_[it->id_ass].eta);
    }
  //getchar();
    }
}

void GenericAnalysis::fill_histos()
{
  // Initi Histos
  TH2F*	hext2=(TH2F*) theRootHandler_->GetTH2("Ext2");
  TH1F* hptgen=(TH1F*)theRootHandler_->GetTH1("PtGen");
  TH1F*  hphigen=(TH1F*)theRootHandler_->GetTH1("PhiGen");
  TH1F* hncout=(TH1F*)theRootHandler_->GetTH1("ngen");
  TH1F*  hnfake=(TH1F*)theRootHandler_->GetTH1("nfake");
  TH1F*  hntracks=(TH1F*)theRootHandler_->GetTH1("ntracks");
  TH1F*  hnpatterns=(TH1F*)theRootHandler_->GetTH1("npatterns");
  TH1F*  hnmatches=(TH1F*)theRootHandler_->GetTH1("nmatches");
  TH1F*  hnmisses=(TH1F*)theRootHandler_->GetTH1("nmisses");

  TH1F* hnstub=(TH1F*)theRootHandler_->GetTH1("nstub");

  TH1F* hpthough=(TH1F*)theRootHandler_->GetTH1("PtHough");
  TH1F* hphihough=(TH1F*)theRootHandler_->GetTH1("PhiHough");
  TH1F* hdpt=(TH1F*)theRootHandler_->GetTH1("deltaPt");
  TH1F*  hdphi=(TH1F*)theRootHandler_->GetTH1("deltaPhi");
  TH2F*	hpt2=(TH2F*) theRootHandler_->GetTH2("Pt_Hough_Gen");
  TH2F*	hphi2=(TH2F*) theRootHandler_->GetTH2("Phi_Hough_Gen");
  TH2F*	hfake2=(TH2F*) theRootHandler_->GetTH2("Fake2");
  TH2F*	hfake2e=(TH2F*) theRootHandler_->GetTH2("Fake2e");
  
  TH2F*	hmiss2=(TH2F*) theRootHandler_->GetTH2("Miss2");
  TH2F*	hfound2=(TH2F*) theRootHandler_->GetTH2("Found2");
  TH2F*	hfound2e=(TH2F*) theRootHandler_->GetTH2("Found2e");
  TH1F*	herrors=(TH1F*)theRootHandler_->GetTH1("Errors");
  TProfile* hpt2p= (TProfile*)theRootHandler_->GetTH1("Pt_Hough_Gen_Profile");
  TProfile* hphi2p=(TProfile*)theRootHandler_->GetTH1("Phi_Hough_Gen_Profile");
  uint32_t nfaketot=0,nfaketot1=0,nfake=0;

  for (std::vector <mctrack_t>::iterator ihbp=theHoughCandidateVector_.begin();ihbp<theHoughCandidateVector_.end();ihbp++)
    {
	
      double pth=(*ihbp).pt;
      double phih=(*ihbp).phi;
      if (pth<1.9) continue;
      hpthough->Fill(log(pth));
      hphihough->Fill(tan(phih));
      double dphmin=9999.;double dptmin=9999.,errmin=9999.;
      std::map<int32_t,mctrack_t>::iterator ismin=theMCMap_.end();
          
		   
					
      for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
	if ( im->second.valid) 
	  {
	    double dph=tan(im->second.phi-phih);						
	    //DEBUG_PRINT(logFile_,"%f %f => %f  \n",tan(im->second.phi),tan(phih),tan(im->second.phi-phih));
	    if (abs((im->second.pt-pth)/pth*dph)>8*3E-4) continue;
	    //if (abs((im->second.pt-pth)/pth)>0.25) continue;
	    // if (abs((im->second.phi-phih))>0.06) continue;
	    //double dph=im->second.phi-phih;
	    double dpt=im->second.pt-pth;
	    //if (abs(dpt)<=dptmin) dptmin=abs(dpt);
	    //if (abs(dph)<=dphmin)	dphmin=abs(dph);

	    //if (abs(dph)>4*5E-3) continue;
	    //if (abs(dpt)/pth>4*6E-2) continue;
				
	    if (abs((im->second.pt-pth)/pth*dph)<errmin)
	      {
		dptmin=abs(dpt);
		dphmin=abs(dph);
		errmin=abs((im->second.pt-pth)/pth*dph);
		ismin=im;
	      }
	
	  }
      //DEBUG_PRINT(logFile_,"Valid Track %d %f %f \n",im->second.id,im->second.pt,im->second.phi);
    
      if (ismin==theMCMap_.end()) 
	{
	  // nfake++;
	  //DEBUG_PRINT(logFile_,"#hits %d Fake Track %f %f et %f %f \n",nbHitPerPattern[k],pth,phih,dptmin,dphmin);
	  // nfaketot1++;
	  
	  //hfake2e->Fill(dptmin/pth,dphmin);
	  double distm=99999.;
	  for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
	    if (im->second.valid || true) 
	      {
		//DEBUG_PRINT(logFile_,"\t %f %f \n",(im->second.pt-pth)/pth,(im->second.phi-phih));
		//double dph=fmod(im->second.phi-phih,PI);
		double dph=tan(im->second.phi-phih);						

		if (abs((im->second.pt-pth)/pth*(dph))<distm) {
		  distm=abs((im->second.pt-pth)/pth*tan(im->second.phi-phih));
		  double dph=tan(im->second.phi-phih);
		  double dpt=im->second.pt-pth;
		  dptmin=abs(dpt);
		  dphmin=abs(dph);
		}
	      }

	  if (distm>8*3E-4) {
	    herrors->Fill(distm);
	    
	    hfake2e->Fill(dptmin/pth,dphmin);
	    hfake2->Fill(pth,tan(phih));
	  }
	  continue;
	}
      //DEBUG_PRINT(logFile_,"Tracks found (%f,%f) nearest is (%f,%f) \n",pth,phih,ismin->second.pt,ismin->second.phi);
      hdpt->Fill((ismin->second.pt-pth)/pth);
      hdphi->Fill(tan(ismin->second.phi-phih));
      hpt2->Fill(ismin->second.pt,pth);
      hphi2->Fill(tan(ismin->second.phi),tan(phih));
      hpt2p->Fill(ismin->second.pt,pth);
      hphi2p->Fill(tan(ismin->second.phi),tan(phih));
      hfound2->Fill(ismin->second.pt,tan(ismin->second.phi));
      hfound2e->Fill((ismin->second.pt-pth)/pth,tan(ismin->second.phi-phih));
      ismin->second.matches++;
    }

  uint32_t nmisses=0;
  for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
    if (im->second.valid)
      {          
	DEBUG_PRINT(logFile_,"Valid Track %d %f %f %f %f matches %d \n",im->second.id,im->second.pt,im->second.phi,im->second.z0,im->second.eta,im->second.matches);
	hncout->Fill(im->second.maxstubs*1.);
	hnmatches->Fill(im->second.matches*1.);
	if (im->second.matches==0)
	  {
	    hmiss2->Fill(im->second.pt,tan(im->second.phi));
	    // for (std::vector<mctrack_t>::iterator ih=theHoughCandidateVector_.begin();ih!=theHoughCandidateVector_.end();ih++)
	    // {
	    // DEBUG_PRINT(logFile_,"Fake %f %f distance %f %f \n",ih->pt,ih->phi,im->second.pt-ih->pt,im->second.phi-ih->phi);
	    // }
	    //getchar();
	    nmisses++;
	  }
				
				
				
      }
  //nfaketot+=nfake;
  DEBUG_PRINT(logFile_,"Number of fake %d  %d %d\n",nfake,nfaketot,nfaketot1); 
  //      getchar();
  if (nfake!=0)hnfake->Fill(nfake*1.);
  hnmisses->Fill(nmisses*1.);

}

void GenericAnalysis::associate()
{
  ngoodmc_=0;
  nmiss_=0;
  nfake_=0;
  
  for (std::vector <mctrack_t>::iterator ihbp=theHoughCandidateVector_.begin();ihbp<theHoughCandidateVector_.end();ihbp++)
    {
	
      double pth=(*ihbp).pt;
      double phih=(*ihbp).phi;
      INFO_PRINT("  Candiddate Hough Track %f %f %d \n",pth,phih,(*ihbp).nhits);
      
      if (pth<thePtCut_*0.9) continue;
      
      double dphmin=9999.;double dptmin=9999.,errmin=9999.;
      std::map<int32_t,mctrack_t>::iterator ismin=theMCMap_.end();
          		   				
      for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
	if (  true ) //im->second.valid ||
	  {
	    double dph=tan(im->second.phi-phih);						
	    double dpt=im->second.pt-pth;
	    if (abs(dph)>1E-2) continue;
	    if (abs(dpt)/pth>1E-1) continue;
	    if (abs(dpt/pth*dph)<errmin)
	      {
		dptmin=abs(dpt);
		dphmin=abs(dph);
		errmin=abs((dpt)/pth*dph);
		ismin=im;
	      }
	
	  }
      if (ismin!=theMCMap_.end())
	{

	  mctrack_t& tk=(*ihbp);
	  //DEBUG_PRINT(logFile_,"==================> Studying %f %f %d \n",tk.pt,tan(tk.phi),theAssociatedTracks_.size());
	  double err_ass=9999.;
	  std::vector<mctrack_t>::iterator isbest=theAssociatedTracks_.end();
	  bool found=false;
	  for (std::vector <mctrack_t>::iterator ia=theAssociatedTracks_.begin();ia<theAssociatedTracks_.end();)
	    {
	      double dph=(ismin->second.phi-(*ia).phi);						
	      double dpt=ismin->second.pt-(*ia).pt;
	      if (abs(dph)>5E-2) {ia++;continue;}
	      if (abs(dpt)/(*ia).pt>1.5E-1) {ia++;continue;}
	      double err=abs(dph);
	      if (err<=errmin)
		{
		  errmin=err;
		  found=true;
		  ia++; continue;
		}
	      else
		{
		  //DEBUG_PRINT(logFile_,"erasing %f %f => %f %f +++ %f %f \n",dph,dpt/(*ia).pt,ismin->second.pt,(*ia).pt,tan(ismin->second.phi),tan((*ia).phi));
		  theAssociatedTracks_.erase(ia);
		}
	    }
	  if (!found)
	    {
	      theAssociatedTracks_.push_back(tk);

	      // DEBUG_PRINT(logFile_,"Adding %g M %g Valid Track %d Pt  %f %f  tphi %f %f \n",err_ass,errmin,ismin->second.id,ismin->second.pt,tk.pt,tan(ismin->second.phi),tan(tk.phi));

	    }
	  ismin->second.matches++;
	  continue;
	}
    
      if (ismin==theMCMap_.end()) 
	{
	  // Not found in MC tracks
	  mctrack_t& tk=(*ihbp);
	  double err_ass=9999.;
	  bool found=false;
	  for (std::vector <mctrack_t>::iterator ia=theFakeTracks_.begin();ia<theFakeTracks_.end();ia++)
	    {
	      double dph=(tk.phi-(*ia).phi);						
	      double dpt=tk.pt-(*ia).pt;
	      if (abs(dph)>5E-2) {ia++;continue;}
	      if (abs(dpt)/(*ia).pt>5E-2) {ia++;continue;}
	      found=true;
	      printf("?????? Bad Associated Hough Track %f %f %d\n",tk.pt,tan(tk.phi),tk.nhits);
	      break;
	    }
	  if (!found && tk.pt>thePtCut_)
	    {
	      theFakeTracks_.push_back(tk);

	      // DEBUG_PRINT(logFile_,"Adding %g M %g Valid Track %d Pt  %f %f  tphi %f %f \n",err_ass,errmin,ismin->second.id,ismin->second.pt,tk.pt,tan(ismin->second.phi),tan(tk.phi));

	    }
	  if (!found && tk.pt<thePtCut_)
	    printf("?????? Bad Hough Track %f %f %d\n",tk.pt,tan(tk.phi),tk.nhits);
	}
      //DEBUG_PRINT(logFile_,"Tracks found (%f,%f) nearest is (%f,%f) \n",pth,phih,ismin->second.pt,ismin->second.phi);
    }

  printf(" Number of Hough candidate %ld and good %ld  and fake %ld \n",theHoughCandidateVector_.size(),theAssociatedTracks_.size(),theFakeTracks_.size());
  for (std::vector <mctrack_t>::iterator ihbp=theAssociatedTracks_.begin();ihbp<theAssociatedTracks_.end();ihbp++)
    printf("Valid Hough Track %f %f %d \n",(*ihbp).pt,tan((*ihbp).phi),(*ihbp).nhits);
  for (std::vector <mctrack_t>::iterator ihbp=theFakeTracks_.begin();ihbp<theFakeTracks_.end();ihbp++)
    printf("-------------> Fake Hough Track %f %f %f %d \n",(*ihbp).pt,tan((*ihbp).phi),(*ihbp).rho0,(*ihbp).nhits);
  nfake_+=theFakeTracks_.size();
  for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
    {
      if (im->second.valid)
	{
	  ngoodmc_++;
	  if (im->second.matches==0)
	    {
	      nmiss_++;
	      printf("++++++ missed MC Track %f %f \n",im->second.pt,tan(im->second.phi));
	    }
	   
	}
    }

  DEBUG_PRINT(logFile_,"Good MC= %d  Missed = %d Fake = %d \n",ngoodmc_,nmiss_,nfake_);
  printf("Good MC= %d  Missed = %d Fake = %d \n",ngoodmc_,nmiss_,nfake_);
  sectmap_[theSector_].goodmc+=ngoodmc_;
  sectmap_[theSector_].missed+=nmiss_;
  sectmap_[theSector_].fake+=nfake_;

  for (uint32_t isect=1;isect<57;isect++)
    if (sectmap_[isect].goodmc)
      printf("%d %d %d %d \n",isect,sectmap_[isect].goodmc,sectmap_[isect].missed,sectmap_[isect].fake);
}

void GenericAnalysis::event_hough()
{
  theHoughCandidateVector_.clear();
  theAssociatedTracks_.clear();
  theFakeTracks_.clear();
  //  HOUGHLOCAL* htl = new HOUGHLOCAL(-PI/2,0.,-0.01,0.01,96,48);
  theHoughLow_->initialise(-PI/2,PI/2,-0.05,0.05,128,theNBinRho_); //160 avant
  theHoughR_->initialise(-PI/2,PI/2,-150.,150.,32,32); //160 avant
  for (std::map<uint32_t,stub_t>::iterator is=theStubMap_.begin();is!=theStubMap_.end();is++)
    {
      //theHoughLow_->fill(is->second.x,is->second.y);
      theHoughLow_->addStub(is->second);
    }
  std::vector < std::pair<double,double> > hbins; hbins.clear();
  std::vector < std::pair<uint16_t,uint16_t> > hibins; hibins.clear();

  uint32_t votemax=(theHoughLow_->getVoteMax()>8)?theHoughLow_->getVoteMax()-4:4;
  
  uint32_t votemin=3;
  //if (theHoughLow_->getVoteMax()>6) votemax=5;
  std::vector<uint32_t> vids;vids.clear();
  theHoughLow_->findMaximumBins(hbins,votemax,&vids);
  //theHoughLow_->findThresholdBins(hbins,9);
  if (theHoughLow_->getVoteMax()<4) return;
  INFO_PRINT(logFile_,"SIZE OF HITS %d\n",(int) vids.size());
  printf("SIZE OF HITS %d %d %d\n",(int) vids.size(),(int) hbins.size(),theHoughLow_->getVoteMax());
#undef  DO_DRAW
#ifdef DO_DRAW

  theHoughLow_->draw(theRootHandler_);
#undef DO_DRAW
#endif
  for (uint32_t i=0;i<theHoughLow_->getNbinTheta();i++)
    for (uint32_t j=0;j<theHoughLow_->getNbinR();j++)
      {
	double R=1./2./TMath::Abs(theHoughLow_->getR(j));
	double pth=0.3*3.8*R/100.;
	if (pth<thePtCut_-0.5) continue;
	std::vector<uint32_t> vid;vid.clear();
	if (theHoughLow_->getHoughImage(i,j)<3) continue;

	std::bitset<24> planes(0);
	if (i>0 && i<theHoughLow_->getNbinTheta()-1 &&j>0 && j<theHoughLow_->getNbinR()-1)
	  {
	    bool nmax=false;
	    for (int ic=-1;ic<=1;ic++)
	      for (int jc=-1;jc<=1;jc++)
		{
		  nmax= (nmax || theHoughLow_->getHoughImage(i,j)<theHoughLow_->getHoughImage(i+ic,j+jc));
		  if (!nmax) 
		    {
#ifndef USE_HV
		      std::vector<uint32_t> v=theHoughLow_->getHoughMap(i+ic,j+jc);
		      for (std::vector<uint32_t>::iterator it=v.begin();it!=v.end();it++)
			if (std::find(vid.begin(),vid.end(), (*it))==vid.end()) {vid.push_back((*it));planes.set(theStubMap_[(*it)].layer,true);}
#else
		      HoughVector v=theHoughLow_->getHoughMap(i+ic,j+jc);
		      for (uint32_t it=0;it<v.size();it++)
			if (std::find(vid.begin(),vid.end(),v[it])==vid.end()) {vid.push_back(v[it]);planes.set(theStubMap_[(*it)].layer,true);}

#endif
		    }
		}
	    if (nmax) continue;
	  

	  }
	else
	  {
#ifndef USE_HV
	    std::vector<uint32_t> v=theHoughLow_->getHoughMap(i,j);
	    for (std::vector<uint32_t>::iterator it=v.begin();it!=v.end();it++)
	      if (std::find(vid.begin(),vid.end(), (*it))==vid.end()) {vid.push_back((*it));planes.set(theStubMap_[(*it)].layer,true);}
#else
	    HoughVector v=theHoughLow_->getHoughMap(i,j);
	    for (uint32_t it=0;it<v.size();it++)
	      if (std::find(vid.begin(),vid.end(),v[it])==vid.end()) {vid.push_back(v[it]);planes.set(theStubMap_[(*it)].layer,true);}

#endif

	  }
	uint32_t nafter=0;
	for (uint32_t ib=0;ib<24;ib++)
	  if (planes[ib]!=0) nafter++;
	if (nafter<3) continue; //was 3
	unsigned long long planhit=planes.to_ulong();
	//printf("%lx \n",planhit);
	if (!planes[5] && endcap) continue;
	if ( barrel && (planhit&0xa0)==0 && (planhit&0xc0)==0 && (planhit&0x60)==0) continue;
	if ( inter && (planhit&0xa0)==0 && (planhit&0xc0)==0 && (planhit&0x60)==0) continue;
	//if (!((planes[5]&&planes[6]) || (planes[7]&&planes[6])||(planes[5]&&planes[7])) )continue;
	
	
	//printf("BIn selected %d,%d  counts %d %f\n",i,j,(int) vid.size(),pth);
#define OLD_VERSION
#ifdef OLD_VERSION
	uint32_t nbinf=64;
	// <5,5-10,10-30,>30
	if (pth<3) nbinf=64;
	if (pth>=3 && pth<5) nbinf=128; // 128
	if (pth>=5  && pth<15) nbinf=128;//192
	if (pth>=15 && pth<=30) nbinf=128;//256
	if (pth>=30 ) nbinf=256;//384

	nbinf /=1; //2 avant
	if (endcap) nbinf/=2;
	uint32_t nbinr=nbinf;
	if (vid.size()>20) nbinf=2*nbinf;

	//if (vid.size()>100) {nbinf=2*nbinf;;nbinr=nbinf;}
#else
	uint32_t nbinf=64;
	// <5,5-10,10-30,>30
	if (pth<3) nbinf=64;
	if (pth>=3 && pth<5) nbinf=128;
	if (pth>=5  && pth<15) nbinf=128;
	if (pth>=15 && pth<=30) nbinf=256;
	if (pth>=30 ) nbinf=512;


#endif
	float ndel=theNDelta_; // 2 avant
	//HOUGHLOCAL *htp = new HOUGHLOCAL(theHoughLow_->getTheta(i)-2*theHoughLow_->getThetaBin(),theHoughLow_->getTheta(i)+2*theHoughLow_->getThetaBin(),theHoughLow_->getR(j)-2*theHoughLow_->getRBin(),theHoughLow_->getR(j)+2*theHoughLow_->getRBin(),nbinf,nbinf);
	theHoughPrecise_->initialise(theHoughLow_->getTheta(i)-ndel*theHoughLow_->getThetaBin(),theHoughLow_->getTheta(i)+ndel*theHoughLow_->getThetaBin(),theHoughLow_->getR(j)-ndel*theHoughLow_->getRBin(),theHoughLow_->getR(j)+ndel*theHoughLow_->getRBin(),nbinf,nbinr);
	theHoughPrecise_->clear();
	printf(" From LowCandidat %f %d Look for bin  val= %x ns %d ith %d ir %d %f %f %f %f %d %d \n",pth,theHoughLow_->getVoteMax(),0xFF,theHoughLow_->getHoughImage(i,j),i,j,theHoughLow_->getTheta(i)-ndel*theHoughLow_->getThetaBin(),theHoughLow_->getTheta(i)+ndel*theHoughLow_->getThetaBin(),theHoughLow_->getR(j)-ndel*theHoughLow_->getRBin(),theHoughLow_->getR(j)+ndel*theHoughLow_->getRBin(),nbinf,nbinr);
	//	printf("LOW %d %d  Size %d PT %f \n",i,j,vid.size(),pth);
	//	getchar();
	for ( std::vector<uint32_t>::iterator itd=vid.begin();itd!=vid.end();itd++)
	  {
	    theHoughPrecise_->addStub(theStubMap_[(*itd)]);
	  }
	analyzePrecise();
	continue;
	//	theHoughPrecise_->draw(d);
	std::vector< std::pair<double,double> > hfbins;hfbins.clear();
	uint32_t vmax=0,nmax=0,imax=0,jmax=0;double thetamax,rmax;
	for (uint32_t ii=0;ii<theHoughPrecise_->getNbinTheta();ii++)
	  for (uint32_t jj=0;jj<theHoughPrecise_->getNbinR();jj++)
	    {
	      
	      if (theHoughPrecise_->getHoughImage(ii,jj)<4) continue;
	      uint neigh=0;
	      double theta=0,r=0;
	     
	      for (int ic=-1;ic<=1;ic++)
		for (int jc=-1;jc<=1;jc++)
		  {
		    if (ic+ii<0) continue;
		    if (ic+ii>=theHoughPrecise_->getNbinTheta()) continue;
		    if (jc+jj<0) continue;
		    if (jc+jj>=theHoughPrecise_->getNbinR()) continue;
		    neigh+=theHoughPrecise_->getHoughImage(ii+ic,jj+jc);
		    r+=(theHoughPrecise_->getR(jj+jc))*theHoughPrecise_->getHoughImage(ii+ic,jj+jc);
		    theta+=(theHoughPrecise_->getTheta(ii+ic))*theHoughPrecise_->getHoughImage(ii+ic,jj+jc);
		  }
	      if (theHoughPrecise_->getHoughImage(ii,jj)>vmax)
		{
		  vmax=theHoughPrecise_->getHoughImage(ii,jj);
		  nmax=neigh;
		  imax=ii;
		  jmax=jj;
		  thetamax=theta/neigh;
		  rmax=r/neigh;
		}
	      else
		if (theHoughPrecise_->getHoughImage(ii,jj)== vmax && neigh>nmax)
		  {
		    vmax=theHoughPrecise_->getHoughImage(ii,jj);
		    nmax=neigh;
		    imax=ii;
		    jmax=jj;
		    thetamax=theta/neigh;
		    rmax=r/neigh;
		  }

	    }
	if (nmax<4) {//delete htp; was 5
	  continue;}
	std::vector<uint32_t> vn;vn.clear();
	//printf("IMAX %d %d \n",imax,jmax);
	for (int ic=-1;ic<=1;ic++)
	  for (int jc=-1;jc<=1;jc++)
	    {
	      if (ic+imax<0) continue;
	      if (ic+imax>=theHoughPrecise_->getNbinTheta()) continue;
	      if (jc+jmax<0) continue;
	      if (jc+jmax>=theHoughPrecise_->getNbinR()) continue;
	      std::vector<uint32_t> v=theHoughPrecise_->getHoughMap(imax+ic,jmax+jc);
	      for (std::vector<uint32_t>::iterator itv=v.begin();itv!=v.end();itv++)
		if (std::find(vn.begin(),vn.end(),(*itv))==vn.end())
		  vn.push_back((*itv));
	    }
	printf("VN SIZE %ld nmax %d \n",vn.size(),nmax);
	std::bitset<24> planez(0);
#ifndef USE_HV
	std::vector<uint32_t> v=theHoughPrecise_->getHoughMap(imax,jmax);
	for (std::vector<uint32_t>::iterator it=vn.begin();it!=vn.end();it++)
	  planez.set(theStubMap_[(*it)].layer,true);
#else
	HoughVector v=theHoughPrecise_->getHoughMap(imax,jmax);
	for (uint32_t it=0;it<v.size();it++)
	  planez.set(theStubMap_[(*it)].layer,true);
#endif

	//std::cout<<planez<<std::endl;
	uint32_t nafte=0;
	for (uint32_t ib=0;ib<24;ib++)
	  if (planez[ib]!=0) nafte++;
	printf("%d %lx \n",nafte,planez.to_ulong());
        if (nafte<5) { continue;} //@@@@ was 4
	if (endcap && planez[5]==0) continue;
	if (barrel && (!((planez[5]&&planez[6]) || (planez[7]&&planez[6])||(planez[5]&&planez[7])) ))continue;
	if (inter && (!((planez[5]&&planez[6]) || (planez[7]&&planez[6])||(planez[5]&&planez[7])) ))continue;

	mctrack_t t;
	
	HOUGHLOCAL::Convert(thetamax,rmax,&t);
	std::vector<uint32_t>::iterator itm=v.begin();
	float x=theStubMap_[(*itm)].x;
	float y=theStubMap_[(*itm)].y;
	//printf("%f %f %f \n",x,y,t.phi);
	if (x>0 && y>0 && t.phi>PI) t.phi-=PI;
	if (x<0 && y>0 && t.phi>PI) t.phi-=PI;
	if (x<0 && y<0 && t.phi<PI) t.phi+=PI;
	if (x>0 && y<0 && t.phi<PI) t.phi+=PI;

	/* if (x>0 && y<0 && t.phi<0) t.phi+=2*PI;
	   if (x>0 && y>0 && t.phi<0) t.phi+=2*PI;
	*/
	t.nhits=nafte;
	theHoughCandidateVector_.push_back(t);
	//@@@hfbins.clear();
        //@@@std::pair<double,double> p(thetamax,rmax);//theHoughPrecise_->getTheta(imax),theHoughPrecise_->getR(jmax));
	//@@@hfbins.push_back(p);



	// if (theHoughPrecise_->getVoteMax()<3) {delete htp;continue;}
	//if (draw) theHoughPrecise_->draw(&rootHandler_);
	//theHoughPrecise_->findMaximumBins(hfbins,3);
        //@@@theHoughPrecise_->draw(theRootHandler_,&hfbins);
	  
	//delete htp;
      }
  //theHoughPrecise_->draw(d,&hfbins);
  /*
    for (std::vector < std::pair<double,double> >::iterator ihb=hbins.begin();ihb<hbins.end();ihb++)
    {
    DEBUG_PRINT(logFile_,"Bins found %f %f => \n",(*ihb).first,(*ihb).second);
    double ith=(*ihb).first;
    double ir=(*ihb).second;
				
    double R=1./2./TMath::Abs(ir);
    double pth=0.3*3.8*R/100.;
    uint32_t nbinf=64;
    if (pth<5) nbinf=128;
    if (pth>=5 && pth<10) nbinf=192;
    if (pth>=10  && pth<30) nbinf=256;
    if (pth>=30) nbinf=320;
    nbinf /=2;
    HOUGHLOCAL *htp = new HOUGHLOCAL(ith-2*htl->getThetaBin(),ith+2*htl->getThetaBin(),ir-2*htl->getRBin(),ir+2*htl->getRBin(),nbinf,nbinf);
				
    htp->clear();
    //cout<<hitIndex<<" Hits found "<<endl;
    //hitf=0;
    //
    //for (std::map<uint32_t,stub_t>::iterator is=theStubMap_.begin();is!=theStubMap_.end();is++)
    htp->fill(is->second.x,is->second.y);
    //
    for ( std::vector<uint32_t>::iterator itd=vids.begin();itd!=vids.end();itd++)
    {
    htp->addStub(theStubMap_[(*itd)]);
    }
	    
    std::vector< std::pair<double,double> > hfbins;hfbins.clear();
    //std::cout<<"OHEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE "<<htp->getVoteMax()<<std::endl;
	    
    if (htp->getVoteMax()<4) {delete htp;continue;}
    //if (draw) htp->draw(&rootHandler_);
    htp->findMaximumBins(hfbins,4);
    htp->draw(d,&hfbins);
	  
    delete htp;
    }
  */

  //delete htl;
}

void GenericAnalysis::analyzePrecise()
{
  std::vector< std::pair<double,double> > hfbins;hfbins.clear();
  uint32_t vmax=0,nmax=0,imax=0,jmax=0;double thetamax,rmax;
  for (uint32_t ii=0;ii<theHoughPrecise_->getNbinTheta();ii++)
    for (uint32_t jj=0;jj<theHoughPrecise_->getNbinR();jj++)
      {
	      
	if (theHoughPrecise_->getHoughImage(ii,jj)<4) continue;
	std::bitset<24> planez(0);
	std::vector<uint32_t> v=theHoughPrecise_->getHoughMap(ii,jj);
	      
	for (std::vector<uint32_t>::iterator it=v.begin();it!=v.end();it++)
	  planez.set(theStubMap_[(*it)].layer,true);
	uint32_t nafte=0;
	for (uint32_t ib=0;ib<24;ib++)
	  if (planez[ib]!=0) nafte++;

        if (nafte<4) { continue;} //@@@@ was 4
	if (!planez[5] && endcap) continue;
	if (barrel && (!((planez[5]&&planez[6]) || (planez[7]&&planez[6])||(planez[5]&&planez[7])) ))continue;
	if (inter && (!((planez[5]&&planez[6]) || (planez[7]&&planez[6])||(planez[5]&&planez[7])) ))continue;
	bool notmax=false;
	for (int ic=-1;ic<=1;ic++)
	  for (int jc=-1;jc<=1;jc++)
	    {
	      if (ic==0 && jc == 0) continue;
	      if (ic+ii<0) continue;
	      if (ic+ii>=theHoughPrecise_->getNbinTheta()) continue;
	      if (jc+jj<0) continue;
	      if (jc+jj>=theHoughPrecise_->getNbinR()) continue;
	      if (theHoughPrecise_->getHoughImage(ii+ic,jj+jc)>theHoughPrecise_->getHoughImage(ii,jj)) 
		{notmax=true;break;}
	    }
	if (notmax) continue;
	

	theHoughR_->clear();
	for (std::vector<uint32_t>::iterator it=v.begin();it!=v.end();it++)
	  theHoughR_->addRStub(theStubMap_[(*it)]);
	//theHoughR_->draw(theRootHandler_);
	if (theHoughR_->getVoteMax()<4) continue;
	double z0=-110.;
	for (uint32_t ir=0;ir<theHoughR_->getNbinTheta();ir++)
	  for (uint32_t jr=0;jr<theHoughR_->getNbinR();jr++)
	    {
	      if (theHoughR_->getHoughImage(ir,jr)==theHoughR_->getVoteMax())
		{
		  double thetar=theHoughR_->getTheta(ir);
		  double rr=theHoughR_->getR(jr);
		  double ar=-1./tan(thetar);
		  double br=rr/sin(thetar);
		  //printf("HOUGH R  %d %f %f %f \n", theHoughR_->getHoughImage(ir,jr),ar,br,-br/ar);
		  z0=-br/ar;
		  //theHoughR_->draw(theRootHandler_);
		  break;
		}
	    }
	uint neigh=0;
	double theta=0,r=0;
	std::vector<uint32_t> vn;vn.clear();
	          
	for (int ic=-1;ic<=1;ic++)
	  for (int jc=-1;jc<=1;jc++)
	    {
	      if (ic+ii<0) continue;
	      if (ic+ii>=theHoughPrecise_->getNbinTheta()) continue;
	      if (jc+jj<0) continue;
	      if (jc+jj>=theHoughPrecise_->getNbinR()) continue;
	      neigh+=theHoughPrecise_->getHoughImage(ii+ic,jj+jc);
	      r+=(theHoughPrecise_->getR(jj+jc))*theHoughPrecise_->getHoughImage(ii+ic,jj+jc);
	      theta+=(theHoughPrecise_->getTheta(ii+ic))*theHoughPrecise_->getHoughImage(ii+ic,jj+jc);


	      std::vector<uint32_t> vp=theHoughPrecise_->getHoughMap(ii+ic,jj+jc);
	      for (std::vector<uint32_t>::iterator itv=vp.begin();itv!=vp.end();itv++)
		if (std::find(vn.begin(),vn.end(),(*itv))==vn.end())
		  vn.push_back((*itv));

	    }
	//printf("%d %d %lx \n",vn.size(),nafte,planez.to_ulong());	
	//if (vn.size()<5) continue;

	theta=theta/neigh;
	r=r/neigh;
	theta=theHoughPrecise_->getTheta(ii);
	r=theHoughPrecise_->getR(jj);



	//std::cout<<planez<<std::endl;
	

	mctrack_t t;
	
	HOUGHLOCAL::Convert(theta,r,&t);
	std::vector<uint32_t>::iterator itm=vn.begin();
	float x=theStubMap_[(*itm)].x;
	float y=theStubMap_[(*itm)].y;
	//printf("%f %f %f \n",x,y,t.phi);
	if (x>0 && y>0 && t.phi>PI) t.phi-=PI;
	if (x<0 && y>0 && t.phi>PI) t.phi-=PI;
	if (x<0 && y<0 && t.phi<PI) t.phi+=PI;
	if (x>0 && y<0 && t.phi<PI) t.phi+=PI;


	t.phierr=theHoughPrecise_->getThetaBin();
	t.pterr=t.pt*theHoughPrecise_->getRBin()/abs(r);
	//printf("%f +- %f  / %f +- %f \n",t.pt,t.pterr,t.phi,t.phierr);
	/* if (x>0 && y<0 && t.phi<0) t.phi+=2*PI;
	   if (x>0 && y>0 && t.phi<0) t.phi+=2*PI;
	*/
	t.nhits=nafte;
	t.z0=z0;
	theHoughCandidateVector_.push_back(t);
		
#ifdef DO_DRAW
	std::vector< std::pair<double,double> > hfbins;hfbins.clear();

        std::pair<double,double> p(theta,r);//theHoughPrecise_->getTheta(imax),theHoughPrecise_->getR(jmax));
	hfbins.push_back(p);



	// if (theHoughPrecise_->getVoteMax()<3) {delete htp;continue;}
	//if (draw) theHoughPrecise_->draw(&rootHandler_);
	//theHoughPrecise_->findMaximumBins(hfbins,3);
        theHoughPrecise_->draw(theRootHandler_,&hfbins);
#endif	  
	//delete htp;
      }

}

void GenericAnalysis::alternativeAssociate()
{
  double PTERRCUT=0.15;//0.15
  double PHIERRCUT=0.05; //5E-2
  ngoodmc_=0;
  nmiss_=0;
  nfake_=0;
  theAssociatedTracks_.clear();
  theFakeTracks_.clear();
  for (std::vector <mctrack_t>::iterator ihbp=theHoughCandidateVector_.begin();ihbp<theHoughCandidateVector_.end();ihbp++)
    {
	
      double pth=(*ihbp).pt;
      double phih=(*ihbp).phi;
      INFO_PRINT("  Candiddate Hough Track %f %f %d \n",pth,phih,(*ihbp).nhits);
      
      //if (pth<thePtCut_*0.9) continue;
      
      double dphmin=9999.;double dptmin=9999.,errmin=9999.;
      std::map<int32_t,mctrack_t>::iterator ismin=theMCMap_.end();
          		   				
      for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
	{
	  double dph=tan(im->second.phi-phih);						
	  double dpt=(im->second.pt-pth);
	  double maxdev=PTERRCUT*pth;
	  if ((*ihbp).pterr>maxdev) maxdev=(*ihbp).pterr;
	  double phidev=PHIERRCUT;
	  if ((*ihbp).phierr>phidev) phidev=(*ihbp).phierr;
	  if (abs(dph)>phidev) continue;
	  
	  if (abs(dpt)>maxdev) continue;
	  //if ((im->second.pt-pth)>(*ihbp).pterr) continue;
	  if (abs(dpt)<errmin)
	    {
	      dptmin=abs(dpt);
	      dphmin=abs(dph);
	      errmin=abs(dpt);
	      ismin=im;
	    }
	
	}
      if (ismin!=theMCMap_.end())
	{

	  mctrack_t& tk=(*ihbp);
	  //DEBUG_PRINT(logFile_,"==================> Studying %f %f %d \n",tk.pt,tan(tk.phi),theAssociatedTracks_.size());
	  double err_ass=9999.;
	  std::vector<mctrack_t>::iterator isbest=theAssociatedTracks_.end();
	  bool found=false;
	  for (std::vector <mctrack_t>::iterator ia=theAssociatedTracks_.begin();ia<theAssociatedTracks_.end();)
	    {
	      double dph=tan(ismin->second.phi-(*ia).phi);						
	      double dpt=(ismin->second.pt-(*ia).pt)/(*ia).pt;

	      double maxdev=PTERRCUT*(*ia).pt;
	      if ((*ia).pterr>maxdev) maxdev=(*ia).pterr;
	      double phidev=PHIERRCUT;
	      if ((*ia).phierr>phidev) phidev=(*ia).phierr;

	      if (abs(dph)>phidev) {ia++;continue;}
	      if (abs(dpt)>maxdev/(*ia).pt) {ia++;continue;}
	      double err=abs(dpt);
	      if (err<=errmin)
		{
		  errmin=err;
		  found=true;
		  ia++; continue;
		}
	      else
		{
		  //DEBUG_PRINT(logFile_,"erasing %f %f => %f %f +++ %f %f \n",dph,dpt/(*ia).pt,ismin->second.pt,(*ia).pt,tan(ismin->second.phi),tan((*ia).phi));
		  theAssociatedTracks_.erase(ia);
		}
	    }
	  if (!found)
	    {
	      tk.id_ass=ismin->first;
	      theAssociatedTracks_.push_back(tk);

	      // DEBUG_PRINT(logFile_,"Adding %g M %g Valid Track %d Pt  %f %f  tphi %f %f \n",err_ass,errmin,ismin->second.id,ismin->second.pt,tk.pt,tan(ismin->second.phi),tan(tk.phi));

	    }
	  ismin->second.matches++;
	  continue;
	}
    
      if (ismin==theMCMap_.end()) 
	{
	  
	  // Not found in MC tracks
	  mctrack_t& tk=(*ihbp);
	  double err_ass=9999.;
	  bool found=false;
	  //  printf("Studying unassociated Hough Track %f %f %d %d\n",tk.pt,tan(tk.phi),tk.nhits,theFakeTracks_.size());
	  for (std::vector <mctrack_t>::iterator ia=theFakeTracks_.begin();ia<theFakeTracks_.end();ia++)
	    {
	      double dph=tan(tk.phi-(*ia).phi);						
	      double dpt=tk.pt-(*ia).pt;
	      double maxdev=PTERRCUT*(*ia).pt;
	      if ((*ia).pterr>maxdev) maxdev=(*ia).pterr;
	      double phidev=PHIERRCUT;
	      if ((*ia).phierr>phidev) phidev=(*ia).phierr;

	      //  printf("%f %f \n",abs(dph),abs(dpt)/(*ia).pt);
	      if (abs(dph)>phidev) {continue;}
	      if (abs(dpt)/(*ia).pt>maxdev/(*ia).pterr) {continue;}
	      found=true;
	      //printf("already associated \n");
	      break;
	    }
	  if (!found && tk.pt>thePtCut_)
	    {
	      //printf("Adding it \n");
	      theFakeTracks_.push_back(tk);

	      // DEBUG_PRINT(logFile_,"Adding %g M %g Valid Track %d Pt  %f %f  tphi %f %f \n",err_ass,errmin,ismin->second.id,ismin->second.pt,tk.pt,tan(ismin->second.phi),tan(tk.phi));

	    }
	  if (!found && tk.pt<thePtCut_)
	    INFO_PRINT("Too low Pt to be added %f %f %d\n",tk.pt,tan(tk.phi),tk.nhits);
	}
      //DEBUG_PRINT(logFile_,"Tracks found (%f,%f) nearest is (%f,%f) \n",pth,phih,ismin->second.pt,ismin->second.phi);
    }

  INFO_PRINT(" Number of Hough candidate %ld and good %ld  and fake %ld \n",theHoughCandidateVector_.size(),theAssociatedTracks_.size(),theFakeTracks_.size());
  for (std::vector <mctrack_t>::iterator ihbp=theAssociatedTracks_.begin();ihbp<theAssociatedTracks_.end();ihbp++)
    INFO_PRINT("Valid Hough Track %f %f %d %f\n",(*ihbp).pt,tan((*ihbp).phi),(*ihbp).nhits,(*ihbp).rho0);
  for (std::vector <mctrack_t>::iterator ihbp=theFakeTracks_.begin();ihbp<theFakeTracks_.end();ihbp++)
    INFO_PRINT("-------------> Fake Hough Track %f %f %f %f %d\n",(*ihbp).pt,tan((*ihbp).phi),(*ihbp).phi,(*ihbp).rho0,(*ihbp).nhits);
  nfake_+=theFakeTracks_.size();
  for (std::map<int32_t,mctrack_t>::iterator im=theMCMap_.begin();im!=theMCMap_.end();im++)
    {
      if (im->second.valid)
	{
	  ngoodmc_++;
	  if (im->second.matches==0)
	    {
	      nmiss_++;
	      INFO_PRINT("++++++ missed MC Track %f %f %f \n",im->second.pt,tan(im->second.phi),im->second.phi);
	    }
	   
	}
    }
  //getchar();
  DEBUG_PRINT(logFile_,"Good MC= %d  Missed = %d Fake = %d \n",ngoodmc_,nmiss_,nfake_);
  INFO_PRINT("Stubs %ld Good MC= %d  Missed = %d Fake = %d \n",theStubMap_.size(),ngoodmc_,nmiss_,nfake_);
  sectmap_[theSector_].goodmc+=ngoodmc_;
  sectmap_[theSector_].missed+=nmiss_;
  sectmap_[theSector_].fake+=nfake_;
}
void GenericAnalysis::PrintSectorMap()
{
  uint32_t ng=0,nm=0,nf=0;
  for (uint32_t isect=1;isect<57;isect++)
    if (sectmap_[isect].goodmc>0)
      {
	float eff=100.*(1.-sectmap_[isect].missed*1./sectmap_[isect].goodmc);
	float fake=100.*(sectmap_[isect].fake*1./sectmap_[isect].goodmc);
	ng+=sectmap_[isect].goodmc;nm+=sectmap_[isect].missed;nf+=sectmap_[isect].fake;
	printf("|%d| %d| %d| %d| %5.2f | % 5.2f| \n",isect,sectmap_[isect].goodmc,sectmap_[isect].missed,sectmap_[isect].fake,eff,fake);
      }
  float eff=100.*(1.-nm*1./ng);
  float fake=100.*(nf*1./ng);
  
  printf("|%d| %d| %d| %d| %5.2f | % 5.2f| \n",255,ng,nm,nf,eff,fake);

}

#ifdef USE_CUDA
void GenericAnalysis::drawph(houghParam* p,DCHistogramHandler* r)
{

  unsigned int h_hough[p->ntheta*p->nrho];
  copyHoughImage(p,h_hough);
  //TH2F* g_hough=(TH2F*) r->GetTH2("GPUHough");
  printf("drawph==> %d %f %f %d %f %f \n",p->ntheta,p->thetamin,p->thetamax,p->nrho,p->rmin,p->rmax);
  TH2F* g_hough=new TH2F("GPUHough","GPUHough",p->ntheta,p->thetamin,p->thetamax,p->nrho,p->rmin,p->rmax);
  for (int ith=0;ith<p->ntheta;ith++)
    for (int ir=0;ir<p->nrho;ir++)
      g_hough->SetBinContent(ith+1,ir+1,h_hough[ith*p->nrho+ir]*1.);
  if (CanvasGPU==NULL)
    {
      CanvasGPU=new TCanvas("CanvasGPU","hough",800,900);
      CanvasGPU->Modified();
      CanvasGPU->Draw();
    }
  CanvasGPU->cd();
  g_hough->Draw("COLZ");
  CanvasGPU->Modified();
  CanvasGPU->Draw();
  //CanvasGPU->WaitPrimitive();
  
  CanvasGPU->Update();
  //CanvasGPU->WaitPrimitive();
  char c;c=getchar();putchar(c); if (c=='.') exit(0);
  delete g_hough;
}
#endif
