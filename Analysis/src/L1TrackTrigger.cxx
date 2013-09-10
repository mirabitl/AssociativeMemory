#define L1TrackTrigger_cxx
#include "L1TrackTrigger.h"
#include "DCHistogramHandler.h"

#include "HoughCartesian.h"

#include "HoughRZ.h"
#include <bitset>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
using namespace std;


std::map<uint32_t,stub_t> stubmap;
std::map<int32_t,mctrack_t> mcmap;
std::vector<pattern_t> patternvec;
std::vector<mctrack_t> houghc;


std::vector<mctrack_t> tk_asso;
std::vector<mctrack_t> tk_fake;

void L1TrackTrigger::fill_histos(DCHistogramHandler* d)
{
  // Initi Histos
  TH2F*	hext2=(TH2F*) d->GetTH2("Ext2");
  TH1F* hptgen=(TH1F*)d->GetTH1("PtGen");
  TH1F*  hphigen=(TH1F*)d->GetTH1("PhiGen");
  TH1F* hncout=(TH1F*)d->GetTH1("ngen");
  TH1F*  hnfake=(TH1F*)d->GetTH1("nfake");
  TH1F*  hntracks=(TH1F*)d->GetTH1("ntracks");
  TH1F*  hnpatterns=(TH1F*)d->GetTH1("npatterns");
  TH1F*  hnmatches=(TH1F*)d->GetTH1("nmatches");
  TH1F*  hnmisses=(TH1F*)d->GetTH1("nmisses");

  TH1F* hnstub=(TH1F*)d->GetTH1("nstub");

  TH1F* hpthough=(TH1F*)d->GetTH1("PtHough");
  TH1F* hphihough=(TH1F*)d->GetTH1("PhiHough");
  TH1F* hdpt=(TH1F*)d->GetTH1("deltaPt");
  TH1F*  hdphi=(TH1F*)d->GetTH1("deltaPhi");
  TH2F*	hpt2=(TH2F*) d->GetTH2("Pt_Hough_Gen");
  TH2F*	hphi2=(TH2F*) d->GetTH2("Phi_Hough_Gen");
  TH2F*	hfake2=(TH2F*) d->GetTH2("Fake2");
  TH2F*	hfake2e=(TH2F*) d->GetTH2("Fake2e");
  
  TH2F*	hmiss2=(TH2F*) d->GetTH2("Miss2");
  TH2F*	hfound2=(TH2F*) d->GetTH2("Found2");
  TH2F*	hfound2e=(TH2F*) d->GetTH2("Found2e");
  TH1F*	herrors=(TH1F*)d->GetTH1("Errors");
  TProfile* hpt2p= (TProfile*)d->GetTH1("Pt_Hough_Gen_Profile");
  TProfile* hphi2p=(TProfile*)d->GetTH1("Phi_Hough_Gen_Profile");
  uint32_t nfaketot=0,nfaketot1=0,nfake=0;

  for (std::vector <mctrack_t>::iterator ihbp=houghc.begin();ihbp<houghc.end();ihbp++)
    {
	
      double pth=(*ihbp).pt;
      double phih=(*ihbp).phi;
      if (pth<1.9) continue;
      hpthough->Fill(log(pth));
      hphihough->Fill(tan(phih));
      double dphmin=9999.;double dptmin=9999.,errmin=9999.;
      std::map<int32_t,mctrack_t>::iterator ismin=mcmap.end();
          
		   
					
      for (std::map<int32_t,mctrack_t>::iterator im=mcmap.begin();im!=mcmap.end();im++)
	if ( im->second.valid) 
	  {
	    double dph=tan(im->second.phi-phih);						
	    //printf("%f %f => %f  \n",tan(im->second.phi),tan(phih),tan(im->second.phi-phih));
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
		    //printf("Valid Track %d %f %f \n",im->second.id,im->second.pt,im->second.phi);
    
      if (ismin==mcmap.end()) 
	{
	  // nfake++;
	  //printf("#hits %d Fake Track %f %f et %f %f \n",nbHitPerPattern[k],pth,phih,dptmin,dphmin);
	  // nfaketot1++;
	  
	  //hfake2e->Fill(dptmin/pth,dphmin);
	  double distm=99999.;
	  for (std::map<int32_t,mctrack_t>::iterator im=mcmap.begin();im!=mcmap.end();im++)
	    if (im->second.valid || true) 
	      {
		//printf("\t %f %f \n",(im->second.pt-pth)/pth,(im->second.phi-phih));
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
		    //printf("Tracks found (%f,%f) nearest is (%f,%f) \n",pth,phih,ismin->second.pt,ismin->second.phi);
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
      for (std::map<int32_t,mctrack_t>::iterator im=mcmap.begin();im!=mcmap.end();im++)
	if (im->second.valid)
          {          
            printf("Valid Track %d %f %f %f %f matches %d \n",im->second.id,im->second.pt,im->second.phi,im->second.z0,im->second.eta,im->second.matches);
	    hncout->Fill(im->second.maxstubs*1.);
	    hnmatches->Fill(im->second.matches*1.);
	    if (im->second.matches==0)
	      {
		hmiss2->Fill(im->second.pt,tan(im->second.phi));
		// for (std::vector<mctrack_t>::iterator ih=houghc.begin();ih!=houghc.end();ih++)
		// {
		// printf("Fake %f %f distance %f %f \n",ih->pt,ih->phi,im->second.pt-ih->pt,im->second.phi-ih->phi);
		// }
		//getchar();
		nmisses++;
	      }
				
				
				
	  }
      //nfaketot+=nfake;
      printf("Number of fake %d  %d %d\n",nfake,nfaketot,nfaketot1); 
      //      getchar();
      if (nfake!=0)hnfake->Fill(nfake*1.);
      hnmisses->Fill(nmisses*1.);

}

void L1TrackTrigger::associate()
{
for (std::vector <mctrack_t>::iterator ihbp=houghc.begin();ihbp<houghc.end();ihbp++)
    {
	
      double pth=(*ihbp).pt;
      double phih=(*ihbp).phi;
      if (pth<1.83) continue;
      
      double dphmin=9999.;double dptmin=9999.,errmin=9999.;
      std::map<int32_t,mctrack_t>::iterator ismin=mcmap.end();
          		   				
      for (std::map<int32_t,mctrack_t>::iterator im=mcmap.begin();im!=mcmap.end();im++)
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
      if (ismin!=mcmap.end())
	{

	  mctrack_t& tk=(*ihbp);
	  //printf("==================> Studying %f %f %d \n",tk.pt,tan(tk.phi),tk_asso.size());
	  double err_ass=9999.;
	  std::vector<mctrack_t>::iterator isbest=tk_asso.end();
	  bool found=false;
	  for (std::vector <mctrack_t>::iterator ia=tk_asso.begin();ia<tk_asso.end();)
	    {
	      double dph=tan(ismin->second.phi-(*ia).phi);						
	      double dpt=ismin->second.pt-(*ia).pt;
	      if (abs(dph)>1E-2) {ia++;continue;}
	      if (abs(dpt)/(*ia).pt>1E-1) {ia++;continue;}
	      double err=abs(dpt/(*ia).pt*dph);
	      if (err<=errmin)
		{
		  errmin=err;
		  found=true;
		  ia++; continue;
		}
	      else
		{
		  //printf("erasing %f %f => %f %f +++ %f %f \n",dph,dpt/(*ia).pt,ismin->second.pt,(*ia).pt,tan(ismin->second.phi),tan((*ia).phi));
		  tk_asso.erase(ia);
		}
	    }
	  if (!found)
	    {
		  tk_asso.push_back(tk);

		  // printf("Adding %g M %g Valid Track %d Pt  %f %f  tphi %f %f \n",err_ass,errmin,ismin->second.id,ismin->second.pt,tk.pt,tan(ismin->second.phi),tan(tk.phi));

	    }
	  ismin->second.matches++;
	  continue;
	}
    
      if (ismin==mcmap.end()) 
	{
	  // Not found in MC tracks
	  mctrack_t& tk=(*ihbp);
	  double err_ass=9999.;
	  bool found=false;
	  for (std::vector <mctrack_t>::iterator ia=tk_fake.begin();ia<tk_fake.end();ia++)
	  {
	      double dph=tan(tk.phi-(*ia).phi);						
	      double dpt=tk.pt-(*ia).pt;
	      if (abs(dph)>5E-3) {ia++;continue;}
	      if (abs(dpt)/(*ia).pt>5E-2) {ia++;continue;}
	      found=true;
	      break;
	  }
	  if (!found && tk.pt>2.)
	    {
		  tk_fake.push_back(tk);

		  // printf("Adding %g M %g Valid Track %d Pt  %f %f  tphi %f %f \n",err_ass,errmin,ismin->second.id,ismin->second.pt,tk.pt,tan(ismin->second.phi),tan(tk.phi));

	    }
	}
		    //printf("Tracks found (%f,%f) nearest is (%f,%f) \n",pth,phih,ismin->second.pt,ismin->second.phi);
    }

 printf(" Number of Hough candidate %d and good %d  and fake %d \n",houghc.size(),tk_asso.size(),tk_fake.size());
 for (std::vector <mctrack_t>::iterator ihbp=tk_asso.begin();ihbp<tk_asso.end();ihbp++)
   printf("Valid Hough Track %f %f \n",(*ihbp).pt,tan((*ihbp).phi));
 for (std::vector <mctrack_t>::iterator ihbp=tk_fake.begin();ihbp<tk_fake.end();ihbp++)
   printf("-------------> Fake Hough Track %f %f \n",(*ihbp).pt,tan((*ihbp).phi));
 nfake_+=tk_fake.size();
 for (std::map<int32_t,mctrack_t>::iterator im=mcmap.begin();im!=mcmap.end();im++)
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

 printf("Good MC= %d  Missed = %d Fake = %d \n",ngoodmc_,nmiss_,nfake_);
}

void L1TrackTrigger::event_hough(DCHistogramHandler* d)
{
  houghc.clear();
  tk_asso.clear();
  tk_fake.clear();
  //  HOUGHLOCAL* htl = new HOUGHLOCAL(-PI/2,0.,-0.01,0.01,96,48);
  htl_->initialise(-PI/2,PI/2,-0.05,0.05,128,160);
  for (std::map<uint32_t,stub_t>::iterator is=stubmap.begin();is!=stubmap.end();is++)
    {
      //htl_->fill(is->second.x,is->second.y);
      htl_->addStub(is->second);
    }
  std::vector < std::pair<double,double> > hbins; hbins.clear();
  std::vector < std::pair<uint16_t,uint16_t> > hibins; hibins.clear();

  uint32_t votemax=htl_->getVoteMax()-1;
  uint32_t votemin=3;
  //if (htl_->getVoteMax()>6) votemax=5;
  std::vector<uint32_t> vids;vids.clear();
  htl_->findMaximumBins(hbins,4,&vids);
  //htl_->findThresholdBins(hbins,9);
  if (htl_->getVoteMax()<4) return;
  printf("SIZE OF HITS %d\n",(int) vids.size());
  //@ htl_->draw(d);

  for (uint32_t i=0;i<htl_->getNbinTheta();i++)
    for (uint32_t j=0;j<htl_->getNbinR();j++)
      {
	double R=1./2./TMath::Abs(htl_->getR(j));
	double pth=0.3*3.8*R/100.;
	if (pth<1.5) continue;
	std::vector<uint32_t> vid;vid.clear();
	if (htl_->getHoughImage(i,j)<4) continue;

	std::bitset<24> planes(0);
	if (i>0 && i<htl_->getNbinTheta()-1 &&j>0 && j<htl_->getNbinR()-1)
	  {
	    bool nmax=false;
	    for (int ic=-1;ic<=1;ic++)
	      for (int jc=-1;jc<=1;jc++)
		{
		  nmax= (nmax || htl_->getHoughImage(i,j)<htl_->getHoughImage(i+ic,j+jc));
		  if (!nmax) 
		    {
#ifndef USE_HV
		      std::vector<uint32_t> v=htl_->getHoughMap(i+ic,j+jc);
		      for (std::vector<uint32_t>::iterator it=v.begin();it!=v.end();it++)
			if (std::find(vid.begin(),vid.end(), (*it))==vid.end()) {vid.push_back((*it));planes.set(LAYER((*it)),true);}
#else
		      HoughVector v=htl_->getHoughMap(i+ic,j+jc);
		      for (uint32_t it=0;it<v.size();it++)
			if (std::find(vid.begin(),vid.end(),v[it])==vid.end()) {vid.push_back(v[it]);planes.set(LAYER(v[it]),true);}

#endif
		    }
		}
	    if (nmax) continue;
	  

	  }
	else
	  {
#ifndef USE_HV
	    std::vector<uint32_t> v=htl_->getHoughMap(i,j);
	    for (std::vector<uint32_t>::iterator it=v.begin();it!=v.end();it++)
	      if (std::find(vid.begin(),vid.end(), (*it))==vid.end()) {vid.push_back((*it));planes.set(LAYER((*it)),true);}
#else
	    HoughVector v=htl_->getHoughMap(i,j);
	    for (uint32_t it=0;it<v.size();it++)
	      if (std::find(vid.begin(),vid.end(),v[it])==vid.end()) {vid.push_back(v[it]);planes.set(LAYER(v[it]),true);}

#endif

	  }
	uint32_t nafter=0;
	for (uint32_t ib=0;ib<24;ib++)
	  if (planes[ib]!=0) nafter++;
	if (nafter<3) continue;
	
	
	//printf("BIn selected %d,%d  counts %d %f\n",i,j,(int) vid.size(),pth);
	uint32_t nbinf=64;
	// <5,5-10,10-30,>30
	if (pth<3) nbinf=64;
	if (pth>=3 && pth<5) nbinf=128;
	if (pth>=5  && pth<15) nbinf=192;
	if (pth>=15 && pth<=30) nbinf=256;
	if (pth>=30 ) nbinf=384;
	nbinf /=2; //2 avant
	//HOUGHLOCAL *htp = new HOUGHLOCAL(htl_->getTheta(i)-2*htl_->getThetaBin(),htl_->getTheta(i)+2*htl_->getThetaBin(),htl_->getR(j)-2*htl_->getRBin(),htl_->getR(j)+2*htl_->getRBin(),nbinf,nbinf);
	htp_->initialise(htl_->getTheta(i)-2*htl_->getThetaBin(),htl_->getTheta(i)+2*htl_->getThetaBin(),htl_->getR(j)-2*htl_->getRBin(),htl_->getR(j)+2*htl_->getRBin(),nbinf,nbinf);
	htp_->clear();
	for ( std::vector<uint32_t>::iterator itd=vid.begin();itd!=vid.end();itd++)
	  {
	    htp_->addStub(stubmap[(*itd)]);
	  }
	//	htp_->draw(d);
	std::vector< std::pair<double,double> > hfbins;hfbins.clear();
	uint32_t vmax=0,nmax=0,imax=0,jmax=0;double thetamax,rmax;
	for (uint32_t ii=0;ii<htp_->getNbinTheta();ii++)
	  for (uint32_t jj=0;jj<htp_->getNbinR();jj++)
	    {
	      
	      if (htp_->getHoughImage(ii,jj)<4) continue;
	      uint neigh=0;
	      double theta=0,r=0;
	      for (int ic=-1;ic<=1;ic++)
		for (int jc=-1;jc<=1;jc++)
		  {
		    if (ic+ii<0) continue;
		    if (ic+ii>=htp_->getNbinTheta()) continue;
		    if (jc+jj<0) continue;
		    if (jc+jj>=htp_->getNbinR()) continue;
		    neigh+=htp_->getHoughImage(ii+ic,jj+jc);
		    r+=(htp_->getR(jj+jc))*htp_->getHoughImage(ii+ic,jj+jc);
		    theta+=(htp_->getTheta(ii+ic))*htp_->getHoughImage(ii+ic,jj+jc);
		  }
	      if (htp_->getHoughImage(ii,jj)>vmax)
		{
		  vmax=htp_->getHoughImage(ii,jj);
		  nmax=neigh;
		  imax=ii;
		  jmax=jj;
		  thetamax=theta/neigh;
		  rmax=r/neigh;
		}
	      else
		if (htp_->getHoughImage(ii,jj)== vmax && neigh>nmax)
		  {
		    vmax=htp_->getHoughImage(ii,jj);
		    nmax=neigh;
		    imax=ii;
		    jmax=jj;
		    thetamax=theta/neigh;
		    rmax=r/neigh;
		  }
	    }
	if (nmax<5) {//delete htp;
	  continue;}
	std::bitset<24> planez(0);
#ifndef USE_HV
	std::vector<uint32_t> v=htp_->getHoughMap(imax,jmax);
	for (std::vector<uint32_t>::iterator it=v.begin();it!=v.end();it++)
	  planez.set(LAYER((*it)),true);
#else
	HoughVector v=htp_->getHoughMap(imax,jmax);
	for (uint32_t it=0;it<v.size();it++)
	  planez.set(LAYER(v[it]),true);
#endif

	//std::cout<<planez<<std::endl;
	uint32_t nafte=0;
	for (uint32_t ib=0;ib<24;ib++)
	  if (planez[ib]!=0) nafte++;
	if (nafte<4) {//delete htp;
	  continue;}
	mctrack_t t;
	HOUGHLOCAL::Convert(thetamax,rmax,&t);
	houghc.push_back(t);
	//@ hfbins.clear();
	//@ std::pair<double,double> p(thetamax,rmax);//htp_->getTheta(imax),htp_->getR(jmax));
	//@ hfbins.push_back(p);
	// if (htp_->getVoteMax()<3) {delete htp;continue;}
	//if (draw) htp_->draw(&rootHandler_);
	//htp_->findMaximumBins(hfbins,3);
	//@ htp_->draw(d,&hfbins);
	  
	//delete htp;
      }
  //htp_->draw(d,&hfbins);
  /*
    for (std::vector < std::pair<double,double> >::iterator ihb=hbins.begin();ihb<hbins.end();ihb++)
    {
    printf("Bins found %f %f => \n",(*ihb).first,(*ihb).second);
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
    //for (std::map<uint32_t,stub_t>::iterator is=stubmap.begin();is!=stubmap.end();is++)
    htp->fill(is->second.x,is->second.y);
    //
    for ( std::vector<uint32_t>::iterator itd=vids.begin();itd!=vids.end();itd++)
    {
    htp->addStub(stubmap[(*itd)]);
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

void L1TrackTrigger::do_ana(std::string fname,uint32_t nevmax)
{
	
  cout<<"On rentre"<<endl;
  cout<<"DCHist created"<<endl;
  htl_ = new HOUGHLOCAL(-PI/2,0.,-0.01,0.01,8,8);
  htp_ = new HOUGHLOCAL(-PI/2,0.,-0.01,0.01,8,8);
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
  TH2F*	hext2=(TH2F*) rootHandler_.BookTH2("Ext2",200,0.,100.,200,-3.,3.);
  TH1F* hptgen=(TH1F*)rootHandler_.BookTH1("PtGen",500,0.5,5.);
  TH1F* hphigen=(TH1F*)rootHandler_.BookTH1("PhiGen",100.,-3.,3.);
  TH1F* hncout=(TH1F*)rootHandler_.BookTH1("ngen",20,0.,20.);
  TH1F*  hnfake=(TH1F*)rootHandler_.BookTH1("nfake",1501,-0.9,1500.1);
  TH1F*  hntracks=(TH1F*)rootHandler_.BookTH1("ntracks",51,-0.9,50.1);
  TH1F*  hnpatterns=(TH1F*)rootHandler_.BookTH1("npatterns",1500,0.,1500.);
  TH1F*  hnmatches=(TH1F*)rootHandler_.BookTH1("nmatches",500,-0.1,499.9);
  TH1F*  hnmisses=(TH1F*)rootHandler_.BookTH1("nmisses",1500,0.,1500.);

  TH1F* hnstub=(TH1F*)rootHandler_.BookTH1("nstub",200,0.,200.);

  TH1F* hpthough=(TH1F*)rootHandler_.BookTH1("PtHough",500,0.5,5.);
  TH1F* hphihough=(TH1F*)rootHandler_.BookTH1("PhiHough",100.,-3.,3.);
  TH1F* hdpt=(TH1F*)rootHandler_.BookTH1("deltaPt",200,-1.,1.);
  TH1F* hdphi=(TH1F*)rootHandler_.BookTH1("deltaPhi",200,-0.05,0.05);
  TH2F*	hpt2=(TH2F*) rootHandler_.BookTH2("Pt_Hough_Gen",200,0.,100.,200,0.,100.);
  TH2F*	hphi2=(TH2F*) rootHandler_.BookTH2("Phi_Hough_Gen",200,-3.,3.,200,-3.,3.);
  TH2F*	hfake2=(TH2F*) rootHandler_.BookTH2("Fake2",200,0.,100.,200,-3.,3.);
  TH2F*	hfake2e=(TH2F*) rootHandler_.BookTH2("Fake2e",200,0.,1.,200,-3.,3.);
  
  TH2F*	hmiss2=(TH2F*) rootHandler_.BookTH2("Miss2",200,0.,100.,200,-3.,3.);
  TH2F*	hfound2=(TH2F*) rootHandler_.BookTH2("Found2",200,0.,100.,200,-3.,3.);
  TH2F*	hfound2e=(TH2F*) rootHandler_.BookTH2("Found2e",200,-0.2,0.2,200,-0.05,0.05);
  TH1F*	herrors=(TH1F*)rootHandler_.BookTH1("Errors",400,0.,8.E-3);	
  TProfile* hpt2p= rootHandler_.BookProfile("Pt_Hough_Gen_Profile",200,0.,100.,0.,100.);
  TProfile* hphi2p=rootHandler_.BookProfile("Phi_Hough_Gen_Profile",200,-3.,3.,-3.,3.);


  /*
    this->initHistograms();

    TH1F* hptgen=(TH1F*)rootHandler_.GetTH1("PtGen");
    TH1F* hphigen=(TH1F*)rootHandler_.GetTH1("PhiGen");
  
    TH1F* hnstub=(TH1F*)rootHandler_.GetTH1("nstub");
    TH1F* hpthough=(TH1F*)rootHandler_.GetTH1("PtHough");
    TH1F* hphihough=(TH1F*)rootHandler_.GetTH1("PhiHough");
    TH1F* hdpt=(TH1F*)rootHandler_.GetTH1("deltaPt");
    TH1F* hdphi=(TH1F*)rootHandler_.GetTH1("deltaPhi");
    TH2F* hpt2=(TH2F*) rootHandler_.GetTH2("Pt_Hough_Gen");
    TH2F* hphi2=(TH2F*) rootHandler_.GetTH2("Phi_Hough_Gen");
    TH2F* hfake2=(TH2F*) rootHandler_.GetTH2("Fake2");
    TH2F* hfake2e=(TH2F*) rootHandler_.GetTH2("Fake2e");
    TH2F* hfound2e=(TH2F*) rootHandler_.GetTH2("Found2e");
    TProfile* hpt2p=(TProfile*) rootHandler_.GetTH1("Pt_Hough_Gen_Profile");
    TProfile* hphi2p=(TProfile*) rootHandler_.GetTH1("Phi_Hough_Gen_Profile");
    TH1F* herrors=(TH1F*)rootHandler_.GetTH1("Errors");

  */
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

  
  if (nevmax!=0 && nevmax<n_entries_MC)
    n_entries_MC=nevmax;
  uint32_t nfaketot=0,nfaketot1=0;
  time_t t0=time(0);
  for (int j=0;j<n_entries_MC;++j){
    PATT->GetEntry(j); // Load entries
    houghc.clear();

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
      mcmap.clear();
      stubmap.clear();
      patternvec.clear();
      TH1F* hnfake=(TH1F*)rootHandler_.GetTH1("nfake");
      TH2F* hmiss2=(TH2F*) rootHandler_.GetTH2("Miss2");
      TH2F* hfound2=(TH2F*) rootHandler_.GetTH2("Found2");
 
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
	    std::map<uint32_t,stub_t>::iterator is=stubmap.find(sid);
	    if (is==stubmap.end())
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
		
		//	printf("%d %d %f \n",stubmap.count(sid),hit_tp[hitIndex],hit_ptGEN[hitIndex]);
		std::pair<uint32_t,stub_t> p(sid,s);
		stubmap.insert(p);
	      }
	    std::map<int32_t,mctrack_t>::iterator im=mcmap.find(hit_tp[hitIndex]);
	    if (mcmap.count(hit_tp[hitIndex])==0 && hit_tp[hitIndex]>=0)
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
		// printf("insert done %d %d \n", mct.id,mcmap.count(mct.id));
		std::pair<uint32_t,mctrack_t> mcp(mct.id,mct);
		mcmap.insert(mcp);

	      }
	    else
	      if (im!=mcmap.end())
		im->second.nstubs|=(1<<hit_layer[hitIndex]);
	    pats.stubs_id.push_back(sid);
	    hitIndex++;

	  }
	
	//printf("Pattern %d Nb stubs %d Total stubs %d Total MC %d \n",k,pats.stubs_id.size(),stubmap.size(),mcmap.size());
	patternvec.push_back(pats);
	//getchar();
      }

      //for (std::map<int32_t,mctrack_t>::iterator im=mcmap.begin();im!=mcmap.end();im++)
      //	im->second.nstubs=0;
      
      for (std::map<int32_t,mctrack_t>::iterator im=mcmap.begin();im!=mcmap.end();im++)
	{


	  uint8_t np=0;
	  for (uint8_t ib=5;ib<=24;ib++)
	    if ((im->second.nstubs>>ib)&1) np++;
	
	  if (im->second.pt<1.2) continue;	  
	  if (np>=4 &&im->second.pt>2 ) 
	    {
	      printf("MC %d NSTUB %x %d PT %f   tg(phi) %f\n",im->second.id,im->second.nstubs,np,im->second.pt,tan(im->second.phi));
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
      for (std::map<int32_t,mctrack_t>::iterator im=mcmap.begin();im!=mcmap.end();im++)
	if (im->second.valid) ngood++;
      printf("MC map size %d Good %d \n",(int) mcmap.size(),ngood);
      event_hough(&rootHandler_);
      associate();
      fill_histos(&rootHandler_);
      continue;
      hntracks->Fill(ngood*1.);
      hitIndex=0;
      htl->clear();
      for(int k=0;k<nb_patterns;k++){

	
	bool draw=true;
	//draw=(j==13);
        htl->clear();
	//htr->clear();
	// affichage du pattern
	// for(int l=0;l<nb_layers;l++){
	// cout<<sstrips[l][k]<<" ";
	// }
	//cout<<" in sector "<<printSector(sector_list[pattern_sector_id[k]])<<endl;
	//cout<<nbHitPerPattern[k]<<" hits in this pattern :"<<k<<endl;
	//affichage des stubs dans le pattern actif
	hitf=hitIndex;
	hitl=hitIndex+nbHitPerPattern[k];
        counts.clear();
	for(int m=0;m<nbHitPerPattern[k];m++)
	  {
	    // cout<<"Layer "<<(int)hit_layer[hitIndex]<<" ladder "<<(int)hit_ladder[hitIndex]<<" Mod "<<(int)hit_zPos[hitIndex]
	    // <<" Segment "<<(int)hit_segment[hitIndex]<<" strip "<<hit_strip[hitIndex]<<" (tp : "<<hit_tp[hitIndex]<<" PT : "<<hit_ptGEN[hitIndex]<<" GeV ETA : "
	    // <<hit_etaGEN[hitIndex]<<" PHI : "<<hit_phiGEN[hitIndex]<<" X:"<<hit_x[hitIndex]<<" Y:"<<hit_y[hitIndex]<<" Z:"<<hit_z[hitIndex]<<")"<<endl;
					
	    htl->fill(hit_x[hitIndex],hit_y[hitIndex]);
	    //htr->fill(sqrt(hit_x[hitIndex]*hit_x[hitIndex]+hit_y[hitIndex]*hit_y[hitIndex]),hit_z[hitIndex]);
	    /*
	      std::map<uint32_t,uint32_t>::iterator itc=counts.find(hit_tp[hitIndex]);
	      if (itc==counts.end())
	      {
              std::pair<uint32_t,uint32_t> p(hit_tp[hitIndex],1);
              counts.insert(p);
	      }
	      else
	      {
              itc->second++;
	      }
	    */
	    hitIndex++;
	  }
	//if (k!=nb_patterns-1) continue;
        uint32_t mcount=0,mindex=0;
        /*
	  for (std::map<uint32_t,uint32_t>::iterator itc=counts.begin();itc!=counts.end();itc++)
          {
	  if (itc->second>mcount)
	  {
	  mcount=itc->second;
	  mindex=itc->first;
	  }
          } 
	  //if (mcount<3) continue;
	  for (uint32_t lidx=hitIndex-1;lidx>0;lidx--)
	  if (hit_tp[lidx]==mindex)
	  {
	  mindex=lidx;break;
	  }
	*/  
	//printf("%d ->%d === %f %f \n",mindex,mcount,hit_ptGEN[mindex],hit_phiGEN[mindex]);
	//if (hit_ptGEN[mindex]>25) continue;
	//if (hit_ptGEN[mindex]<2) continue;
	
	if (hnstub==NULL)
	  {
	   
	   
	  }
	
	//hncout->Fill(mcount*1.);
	hnstub->Fill(nbHitPerPattern[k]);
	//getchar();
	// for(int k=0;k<nb_tracks;k++){
	// cout<<"Track "<<k<<" : PT="<<track_pt[k]<<" PHI0="<<track_phi[k]<<" ETA="<<track_eta[k]<<" D0="<<track_d0[k]<<" Z0="<<track_z0[k]<<endl;
	// }
	//cout<<endl;
	uint32_t votemin=3;
	if (nbHitPerPattern[k]>7) votemin=4;
	if (nbHitPerPattern[k]>8) votemin=5;
	std::vector < std::pair<double,double> > hbins; hbins.clear();
	uint32_t votemax=htl->getVoteMax()-1;
	//if (htl->getVoteMax()>6) votemax=5;
	htl->findMaximumBins(hbins,TMath::Max(votemin,votemax));
	//htl->findThresholdBins(hbins,9);
	if (htl->getVoteMax()<3) continue;
	for (std::vector < std::pair<double,double> >::iterator ihb=hbins.begin();ihb<hbins.end();ihb++)
	  {
	    // printf("Bins found %f %f => \n",(*ihb).first,(*ihb).second);
	  }
	if (draw) htl->draw(&rootHandler_);
	// std::vector < std::pair<double,double> > hrbins; hrbins.clear();
	// htr->findMaximumBins(hrbins,3);//TMath::Max((uint32_t)3, (htr->getVoteMax()-1)));	
	// htr->draw(&rootHandler_,&hrbins);
	//if (nbHitPerPattern[k]>0) continue;

	for (std::vector < std::pair<double,double> >::iterator ihb=hbins.begin();ihb<hbins.end();ihb++)
	  {
	    double ith=(*ihb).first;
	    double ir=(*ihb).second;
				
				
				
	    //printf("Bin  studied %f %f %f %f => \n",ith-2*htl->getThetaBin(),ith+2*htl->getThetaBin(),ir-2*htl->getRBin(),ir+2*htl->getRBin());
	    //HOUGHLOCAL::PrintConvert(ith,ir);
	    double R=1./2./TMath::Abs(ir);
	    double pth=0.3*3.8*R/100.;
	    uint32_t nbinf=64;
	    if (pth<5) nbinf=128;
	    if (pth>=5 && pth<10) nbinf=192;
	    if (pth>=10  && pth<30) nbinf=256;
	    if (pth>=30) nbinf=320;

	    HOUGHLOCAL *htp = new HOUGHLOCAL(ith-2*htl->getThetaBin(),ith+2*htl->getThetaBin(),ir-2*htl->getRBin(),ir+2*htl->getRBin(),nbinf,nbinf);
				
	    htp->clear();
	    //cout<<hitIndex<<" Hits found "<<endl;
	    //hitf=0;
	    for (int is=hitf;is<hitl;is++)
	      htp->fill(hit_x[is],hit_y[is]);

	    std::vector< std::pair<double,double> > hfbins;hfbins.clear();
	    //std::cout<<"OHEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE "<<htp->getVoteMax()<<std::endl;
	    
	    if (htp->getVoteMax()<votemin) {delete htp;continue;}
	    //if (draw) htp->draw(&rootHandler_);
	    htp->findMaximumBins(hfbins,TMath::Max(htp->getVoteMax()-1,votemin));
	    if (draw) htp->draw(&rootHandler_,&hfbins);
	    if(hfbins.size()>0)
	      {
		for (std::vector < std::pair<double,double> >::iterator ihbp=hfbins.begin();ihbp<hfbins.end();ihbp++)
		  {
		    //printf(" Track ===> Theta  %f r %f \n",(*ihbp).first,(*ihbp).second);
		    double theta=(*ihbp).first;
		    double r=(*ihbp).second;
		
		    double a=-1./tan(theta);
		    double b=r/sin(theta);
		
		
		    //
		    double R=1./2./TMath::Abs(r);
		    double pth=0.3*3.8*R/100.;
		    double phih=atan(a);
		    if (pth<2.) continue;
		    hpthough->Fill(pth);
		    hphihough->Fill(a);
		    double dphmin=9999.;double dptmin=9999.,errmin=9999.;
		    std::map<int32_t,mctrack_t>::iterator ismin=mcmap.end();
          
		   
					
		    for (std::map<int32_t,mctrack_t>::iterator im=mcmap.begin();im!=mcmap.end();im++)
		      if ( im->second.valid) 
			{
							

			  if (abs((im->second.pt-pth)/pth*(im->second.phi-phih))>1*3E-4) continue;
			  //if (abs((im->second.pt-pth)/pth)>0.25) continue;
			  // if (abs((im->second.phi-phih))>0.06) continue;
			  double dph=im->second.phi-phih;
			  double dpt=im->second.pt-pth;
			  //if (abs(dpt)<=dptmin) dptmin=abs(dpt);
			  //if (abs(dph)<=dphmin)	dphmin=abs(dph);

			  //if (abs(dph)>4*5E-3) continue;
			  //if (abs(dpt)/pth>4*6E-2) continue;
				
			  if (abs((im->second.pt-pth)/pth*(im->second.phi-phih))<errmin)
			    {
			      dptmin=abs(dpt);
			      dphmin=abs(dph);
			      errmin=abs((im->second.pt-pth)/pth*(im->second.phi-phih));
			      ismin=im;
			    }
	
			}
		    //printf("Valid Track %d %f %f \n",im->second.id,im->second.pt,im->second.phi);
    
		    if (ismin==mcmap.end()) {nfake++;
		      //printf("#hits %d Fake Track %f %f et %f %f \n",nbHitPerPattern[k],pth,phih,dptmin,dphmin);
		      nfaketot1++;
		      hfake2->Fill(pth,phih);
		      //hfake2e->Fill(dptmin/pth,dphmin);
		      double distm=99999.;
		      for (std::map<int32_t,mctrack_t>::iterator im=mcmap.begin();im!=mcmap.end();im++)
			if (im->second.valid) 
			  {
			    //printf("\t %f %f \n",(im->second.pt-pth)/pth,(im->second.phi-phih));
			    if (abs((im->second.pt-pth)/pth*(im->second.phi-phih))<distm) distm=abs((im->second.pt-pth)/pth*(im->second.phi-phih));
			  }
		      herrors->Fill(distm);

		      mctrack_t t;
		      t.pt=pth;
		      t.phi=phih;
		      houghc.push_back(t);
		      //getchar();
		      continue;}
		    //printf("Tracks found (%f,%f) nearest is (%f,%f) \n",pth,phih,ismin->second.pt,ismin->second.phi);
		    hdpt->Fill((ismin->second.pt-pth)/pth);
		    hdphi->Fill(ismin->second.phi-phih);
		    hpt2->Fill(ismin->second.pt,pth);
		    hphi2->Fill(ismin->second.phi,phih);
		    hpt2p->Fill(ismin->second.pt,pth);
		    hphi2p->Fill(ismin->second.phi,phih);
		    hfound2->Fill(ismin->second.pt,ismin->second.phi);
		    hfound2e->Fill((ismin->second.pt-pth)/pth,ismin->second.phi-phih);
		    ismin->second.matches++;
		  }
		if (draw) htp->draw(&rootHandler_,&hfbins);
	      }
	    delete htp;
	  }
	/*
				
	 */
			
      }
      uint32_t nmisses=0;
      for (std::map<int32_t,mctrack_t>::iterator im=mcmap.begin();im!=mcmap.end();im++)
	if (im->second.valid)
          {          
            printf("Valid Track %d %f %f %f %f matches %d \n",im->second.id,im->second.pt,im->second.phi,im->second.z0,im->second.eta,im->second.matches);
	    hncout->Fill(im->second.maxstubs*1.);
	    hnmatches->Fill(im->second.matches*1.);
	    if (im->second.matches==0)
	      {
		hmiss2->Fill(im->second.pt,im->second.phi);
		// for (std::vector<mctrack_t>::iterator ih=houghc.begin();ih!=houghc.end();ih++)
		// {
		// printf("Fake %f %f distance %f %f \n",ih->pt,ih->phi,im->second.pt-ih->pt,im->second.phi-ih->phi);
		// }
		//getchar();
		nmisses++;
	      }
				
				
				
	  }
      nfaketot+=nfake;
      printf("Number of fake %d  %d %d\n",nfake,nfaketot,nfaketot1); 
      //      getchar();
      if (nfake!=0)hnfake->Fill(nfake*1.);
      hnmisses->Fill(nmisses*1.);
      //getchar();
      if (nfake>40)
	{
	  printf("LARGE %d \n",j);
	  //getchar();
	}
    }
		
  }
  delete PATT;
  rootHandler_.writeHistograms("toto.root");
	
}

void L1TrackTrigger::Loop(TApplication* theApp)
{
  //   In a ROOT session, you can do:
  //      Root > .L L1TrackTrigger.C
  //      Root > L1TrackTrigger t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;
	
  Long64_t nentries = fChain->GetEntriesFast();
  //DCHistogramHandler rootHandler_;
  //HOUGHLOCAL* htl = new HOUGHLOCAL(PI/2.,PI,-0.02,0.02,32,32);
  HOUGHLOCAL* htl = new HOUGHLOCAL(PI/2,PI,-0.05,0.05,512,512);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //if (STUB_tp->size()==0) continue;
    htl->clear();
    double Ptmin=99999.;
    double Ptmax=0;

    for (int is=0;is<STUB_n;is++)
      {	
	if (TMath::Abs((*STUB_pt)[is])<1.5) continue;
	if (TMath::Abs((*STUB_deltas)[is])>2.) continue;
	double pt=sqrt((*STUB_pxGEN)[is]*(*STUB_pxGEN)[is]+(*STUB_pyGEN)[is]*(*STUB_pyGEN)[is]);
	if (pt<Ptmin && pt>1.5) Ptmin=pt;
	if (pt>Ptmax && pt>1.5) Ptmin=pt;
	if ((*STUB_x)[is]<0) continue;
	if ((*STUB_y)[is]<0) continue;
	printf("%f %f %f \n",pt,(*STUB_pt)[is],(*STUB_deltas)[is]);
			
	htl->fill((*STUB_x)[is],(*STUB_y)[is]);
      }
    ///if (Ptmin<2 || Ptmax>25) continue;
		

#ifdef OLDWAY
    std::vector < std::pair<uint32_t,uint32_t> > hbins; hbins.clear();
    htl->findMaxima(hbins,(uint32_t) (htl->getVoteMax()-1));
    for (std::vector < std::pair<uint32_t,uint32_t> >::iterator ihb=hbins.begin();ihb<hbins.end();ihb++)
      {
	printf("Bin %d %d => %d \n",(*ihb).first,(*ihb).second,htl->getValue((*ihb).first,(*ihb).second));
      }
    htl->draw(&rootHandler_);
    for (std::vector < std::pair<uint32_t,uint32_t> >::iterator ihb=hbins.begin();ihb<hbins.end();ihb++)
      {
	uint32_t ith=(*ihb).first;
	uint32_t ir=(*ihb).second;
	printf("Bin  studied %d %d => %d  new bin %f %f %f %f \n",ith,ir,htl->getValue(ith,ir),htl->getTheta(ith-1),htl->getTheta(ith+2),htl->getR(ir-1),htl->getR(ir+2));
	HOUGHLOCAL htp(htl->getTheta(ith-1),htl->getTheta(ith+2),htl->getR(ir-1),htl->getR(ir+2),64,64);
	htp.clear();
	//htp.fill(&allpoints_);
	for (int is=0;is<STUB_n;is++)
	  htp.fill((*STUB_x)[is],(*STUB_y)[is]);
	//std::vector < std::pair<uint32_t,uint32_t> > hpbins; hpbins.clear();
	//htp.findMaxima(hpbins,1);
	//for (std::vector < std::pair<uint32_t,uint32_t> >::iterator ihbp=hpbins.begin();ihbp<hpbins.end();ihbp++)
	//{
	//printf("Bin %d %d => %d \n",(*ihbp).first,(*ihbp).second,htp.getValue((*ihbp).first,(*ihbp).second));
	//	}
	//htp.draw(&rootHandler_,&hpbins);
	std::vector< std::pair<double,double> > hfbins;hfbins.clear();
	htp.findMaxima(hfbins,TMath::Max(htp.getVoteMax(),(uint32_t)5));
			
	for (std::vector < std::pair<double,double> >::iterator ihbp=hfbins.begin();ihbp<hfbins.end();ihbp++)
	  {
	    printf(" Track ===> Theta  %f r %f \n",(*ihbp).first,(*ihbp).second);
	  }
	htp.draw(&rootHandler_,&hfbins);
      }
		
#endif

    std::vector < std::pair<double,double> > hbins; hbins.clear();
    printf("Maximal vote %d \n",htl->getVoteMax());
    htl->findMaximumBins(hbins,TMath::Max((uint32_t)5, (htl->getVoteMax()-3)));
    for (std::vector < std::pair<double,double> >::iterator ihb=hbins.begin();ihb<hbins.end();ihb++)
      {
	printf("Bin %f %f => \n",(*ihb).first,(*ihb).second);
      }
    htl->draw(&rootHandler_);
			
    for (std::vector < std::pair<double,double> >::iterator ihb=hbins.begin();ihb<hbins.end();ihb++)
      {
	double ith=(*ihb).first;
	double ir=(*ihb).second;
				
				
				
	//printf("Bin  studied %f %f => \n",ith,ir);
	HOUGHLOCAL::PrintConvert(ith,ir);
			
	HOUGHLOCAL *htp = new HOUGHLOCAL(ith-htl->getThetaBin(),ith+htl->getThetaBin(),ir-htl->getRBin(),ir+htl->getRBin(),32,32);
				
	htp->clear();
	for (int is=0;is<STUB_n;is++)
	  {	
	    if (TMath::Abs((*STUB_pt)[is])<1.5) continue;
	    if (TMath::Abs((*STUB_deltas)[is])>2.) continue;
	    //if ((*STUB_tp)[is]!=0) continue;
	    if ((*STUB_x)[is]<0) continue;
	    if ((*STUB_y)[is]<0) continue;

	    htp->fill((*STUB_x)[is],(*STUB_y)[is]);
	  }	
	std::vector< std::pair<double,double> > hfbins;hfbins.clear();
	htp->findMaximumBins(hfbins,TMath::Max(htp->getVoteMax()-1,(uint32_t)5));
	if(hfbins.size()>0)
	  {
	    for (std::vector < std::pair<double,double> >::iterator ihbp=hfbins.begin();ihbp<hfbins.end();ihbp++)
	      {
		printf(" Track ===> Theta  %f r %f \n",(*ihbp).first,(*ihbp).second);
	      }
	    htp->draw(&rootHandler_,&hfbins);
	  }
	delete htp;
      }
    /*
				
     */
			
			
		
  }
}
