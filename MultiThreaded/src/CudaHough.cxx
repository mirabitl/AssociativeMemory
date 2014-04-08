#include <stdint.h>
#include <stdio.h> 
#include "CudaHough.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <string.h>
#include <TMath.h>
#include "DCHistogramHandler.h"

CudaHough::CudaHough(HoughCut* cuts) :theCuts_(cuts)
{
  theNStub_=0;
  theX_=NULL;
  theY_=NULL;
  theZ_=NULL;
  theLayer_=NULL;

  createHough(&ph_,768,3072,768);
  for (int i=0;i<96;i++)
    createHough(&phcand_[i]);
  for (int i=0;i<64;i++)
    createHough(&phrcand_[i]);
  createStreams(96);
  initialiseTimer();

  createTkletEvent(&et_);
  printf("%d size \n",sizeof(ctklevent));
}
void CudaHough::DefaultCuts()
{
  theCuts_->RhoMin=-0.0031;
  theCuts_->RhoMax=0.0031;
  theCuts_->NRho=6;
  theCuts_->NTheta=48;
  theCuts_->NStubLow=4;
  theCuts_->NLayerRow=4;
  theCuts_->NStubLowCandidate=5;
  theCuts_->NBins3GeV=56;
  theCuts_->NBins5GeV=128; 
  theCuts_->NBins15GeV=128;
  theCuts_->NBins30GeV=128;
  theCuts_->NBins100GeV=128;
  theCuts_->NDelBarrel=1.5;
  theCuts_->NDelInter=2.1;
  theCuts_->NDelEndcap=2.7;
  theCuts_->NStubHigh=5;
  theCuts_->NLayerHigh=5;
  theCuts_->NStubHighCandidate=5;
}
void CudaHough::Compute(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer)
{
  theNStub_=nstub;
  theX_=x;
  theY_=y;
  theZ_=z;  
  theLayer_=layer;
  theCandidateVector_.clear();
  // Initialisation depending on sector 
  bool barrel=isel>=16 && isel<40;
  bool inter=(isel>=8 &&isel<16)||(isel>=40&&isel<48);
  bool endcap=(isel<8)||(isel>=48);
  float thmin=-PI/2,thmax=PI/2;
  float rhmin=theCuts_->RhoMin,rhmax=theCuts_->RhoMax;
  //printf("On appelle le GPU %d \n",theNStub_);
  int ntheta=160;
  int nrho=theCuts_->NRho;//8//12//192;
  //initialiseHough(&ph,gpu_nstub,ntheta,nrho,-PI/2,PI/2,-0.06,0.06);
  if (barrel || inter || endcap)
    {
      ntheta=theCuts_->NTheta;//64;
      
      if (isel%4==0) thmin=1.32;
      if (isel%4==1) thmin=-1.04;
      if (isel%4==2) thmin=-0.24;
      if (isel%4==3) thmin=0.51;
      thmax=thmin+1.25;
    }

  if (theNStub_>400)
    {
      ntheta*=2;
      nrho*=2;
    }
  
  initialiseHough(&ph_,theNStub_,ntheta,nrho,thmin,thmax,rhmin,rhmax);
  // Rough process
  fillConformalHough(&ph_,theX_,theY_,theZ_);
  fillLayerHough(&ph_,theLayer_);
		  //clearHough(&ph);
  //  std::cout<<endcap<<std::endl;
  processHough(&ph_,theCuts_->NStubLow,theCuts_->NLayerRow,0,-1,endcap);
  //printf("SECTOR %d gives %d candidates Max val %d STubs %d\n",isel,ph_.h_cand[0],ph_.max_val,ph_.nstub);
  // Precise HT filling
  uint32_t nc=(int)ph_.h_cand[0];
  if (nc>96) nc=96;
  for (int ic=0;ic<nc;ic++)
    {
      clearHough(&phcand_[ic]);
    }

  // Loop on candidate
  for (int ic=0;ic<nc;ic++)
    {
      phcand_[ic].h_reg[20]=0;
      int pattern=ph_.h_cand[ic+1]; // vcand[ic]
      int ith=pattern&0X3FF;
      int ir=(pattern>>10)&0x3FF;
      //ith=(vcand[ic])&0x3FF;
      //ir=(vcand[ic]>>10)&0x3FF;
      int ns=(pattern>>20)&0x3FF;
      if (ns<theCuts_->NStubLowCandidate) continue;//if (ns<3) continue;
      double PT=1./2./fabs(GET_R_VALUE(ph_,ir))*0.3*3.8/100.;
      if (PT<1.5) continue;
      //printf("%f \n",fabs(GET_R_VALUE(ph,ir)));
      uint32_t nbinf=64;
      // <5,5-10,10-30,>30
      if (PT<3) nbinf=theCuts_->NBins3GeV;
      if (PT>=3 && PT<5) nbinf=theCuts_->NBins5GeV; // 128
      if (PT>=5  && PT<15) nbinf=theCuts_->NBins15GeV;//192
      if (PT>=15 && PT<=30) nbinf=theCuts_->NBins30GeV;//256
      if (PT>=30 ) nbinf=theCuts_->NBins100GeV;//256



      uint32_t nbinr=nbinf;

      if (ns>20 ) nbinf=2*nbinf;
		  
      float ndel=theCuts_->NDelBarrel;
      if (inter) 
	ndel=theCuts_->NDelInter;
      else
	if (endcap)
	  ndel=theCuts_->NDelEndcap;


        


      float tmi=GET_THETA_VALUE(ph_,ith)-ndel*ph_.thetabin;
      
      float tma=GET_THETA_VALUE(ph_,ith)+ndel*ph_.thetabin;
      float rmi=GET_R_VALUE(ph_,ir)-ndel*ph_.rbin;
      float rma=GET_R_VALUE(ph_,ir)+ndel*ph_.rbin;
      
      initialiseHough(&phcand_[ic],theNStub_,nbinf,nbinr,tmi,tma,rmi,rma);	    
      //      std::cout<<endcap<<std::endl;
      copyPositionHough(&ph_,pattern,&phcand_[ic],0,false,ic,endcap);
    }

  synchronize();

		  

		
  //Precise HT processing
		 
  for (int ic=0;ic<nc;ic++)
    {
      if (phcand_[ic].h_reg[20]>0)
	{
	  phcand_[ic].nstub=int( phcand_[ic].h_reg[20]);
	  //	  std::cout<<endcap<<std::endl;
	  processHough(&phcand_[ic],theCuts_->NStubHigh,theCuts_->NLayerHigh,0,ic,endcap);
	  
	}
    }

  synchronize();
  // Finael analysis of High precision candidate

  for (int ic=0;ic<nc;ic++)
    {
      if (phcand_[ic].h_reg[20]>0)
	{
	  uint32_t nch=(int)phcand_[ic].h_cand[0];
	  if (nch>64) nch=64;
	  for (int ici=0;ici<nch;ici++)
	    {
	      int patterni=phcand_[ic].h_cand[ici+1]; 
	      int ithi=patterni&0X3FF;
	      int iri=(patterni>>10)&0x3FF;
	      
	      if (((patterni>>20)&0x3FF)<theCuts_->NStubHighCandidate) continue;

	      mctrack_t t;
	      // RZ  & R Phi regression
	      initialiseHough(&phrcand_[ici],theNStub_,32,32,-PI/2,PI/2,-150.,150.);
	      //	      std::cout<<endcap<<std::endl;
	      copyPositionHough(&phcand_[ic],patterni,&phrcand_[ici],1,true,-1,endcap);
	      phrcand_[ici].nstub=int( phrcand_[ici].h_reg[20]);
	      if (phrcand_[ici].h_reg[60+6]<1.7) continue;
	      if ( phrcand_[ici].h_reg[20]<=0) continue;
			      
	      if ( phrcand_[ici].h_reg[70+9]<1.5) continue; //at least 2 Z points
	      t.z0=-phrcand_[ici].h_reg[70+1]/phrcand_[ici].h_reg[70+0];
	      t.eta=phrcand_[ici].h_reg[70+8];
	      if ( fabs(t.z0)>30.) continue; //30 avant
			  
	      
	      float theta=GET_THETA_VALUE(phcand_[ic],ithi);
	      float r=GET_R_VALUE(phcand_[ic],iri);

	      double a=-1./tan(theta);
	      double b=r/sin(theta);
			      
		
			  //
	      double R=1./2./fabs(r);
	      double xi=-a/2./b;
	      double yi=1./2./b;
	      double g_pt=0.3*3.8*R/100.;
			  //g_phi=atan(a);
	      double g_phi=theta-PI/2.;
	      if (g_phi<0) g_phi+=2*PI;
	      CudaHough::Convert(theta,r,&t);
	      t.nhits=(patterni>>20)&0x3FF;
	      t.theta=theta;
	      t.r=r;

	      t.pt=phrcand_[ici].h_reg[60+6];
	      t.phi=phrcand_[ici].h_reg[60+2];
	      t.nhits=(patterni>>20)&0x3FF;
	      //	      std::cout<<endcap<<std::endl;
	      getChi2(&phrcand_[ici],endcap);
	      t.chi2=phrcand_[ici].h_reg[80];
	      t.chi2z=phrcand_[ici].h_reg[81];


	      t.layers.clear();
	      int layers[2048];
	      copyLayers(&phrcand_[ici],layers);
	      for (int ist=0;ist<phrcand_[ici].nstub;ist++)
		t.layers.push_back(layers[ist]);

	      theCandidateVector_.push_back(t);

			      
	    }
	}
    }
		 









		  
		
		  
  //  printf("Fin du GPU %ld \n",	theCandidateVector_.size() );



}

void CudaHough::ComputeOneShot(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer)
{
  theNStub_=nstub;
  theX_=x;
  theY_=y;
  theZ_=z;  
  theLayer_=layer;
  theCandidateVector_.clear();
  // Initialisation depending on sector 
  bool barrel=isel>=16 && isel<40;
  bool inter=(isel>=8 &&isel<16)||(isel>=40&&isel<48);
  bool endcap=(isel<8)||(isel>=48);
  float thmin=-PI/2,thmax=PI/2;
  float rhmin=theCuts_->RhoMin,rhmax=theCuts_->RhoMax;
  //printf("On appelle le GPU %d \n",theNStub_);
  int ntheta=160;
  int nrho=theCuts_->NRho;//8//12//192;
  //initialiseHough(&ph,gpu_nstub,ntheta,nrho,-PI/2,PI/2,-0.06,0.06);
  if (barrel || inter || endcap )
    {
      ntheta=theCuts_->NTheta;//64;
      if (isel%4==0) thmin=1.32;
      if (isel%4==1) thmin=-1.04;
      if (isel%4==2) thmin=-0.24;
      if (isel%4==3) thmin=0.51;
      thmax=thmin+1.25;
    }

  
  if (theNStub_>400)
    {
      ntheta*=2;
      nrho*=2;
    }
  ntheta=960;
  nrho=156;
  if (inter)
    {
      ntheta=1056;
      nrho=88;
    }
  if (endcap)
    {
      ntheta=1056;
      nrho=64;
    }
  theCuts_->NLayerRow=5;

  initialiseHough(&ph_,theNStub_,ntheta,nrho,thmin,thmax,rhmin,rhmax);
  // Rough process
  fillConformalHough(&ph_,theX_,theY_,theZ_);
  fillLayerHough(&ph_,theLayer_);
		  //clearHough(&ph);
  processHough(&ph_,theCuts_->NStubLow,theCuts_->NLayerRow,0,-1,endcap);
  //printf("SECTOR %d gives %d candidates Max val %d STubs %d\n",isel,ph_.h_cand[0],ph_.max_val,ph_.nstub);
  // Precise HT filling
  uint32_t nc=(int)ph_.h_cand[0];
  if (nc>512) nc=512;
  clearHough(&phcand_[0]);
  

  // Loop on candidate
  for (int ic=0;ic<nc;ic++)
    {
      phcand_[0].h_reg[20]=0;
      int pattern=ph_.h_cand[ic+1]; // vcand[ic]
      int ith=pattern&0X3FF;
      int ir=(pattern>>10)&0x3FF;
      //ith=(vcand[ic])&0x3FF;
      //ir=(vcand[ic]>>10)&0x3FF;
      int ns=(pattern>>20)&0x3FF;
      if (ns<theCuts_->NStubLowCandidate) continue;//if (ns<3) continue;
      double PT=1./2./fabs(GET_R_VALUE(ph_,ir))*0.3*3.8/100.;
      if (PT<1.5) continue;
      //printf("%f \n",fabs(GET_R_VALUE(ph,ir)));

      mctrack_t t;
      // RZ  & R Phi regression
      initialiseHough(&phcand_[0],theNStub_,32,32,-PI/2,PI/2,-150.,150.);
      copyPositionHough(&ph_,pattern,&phcand_[0],1,true,-1,endcap);
      phcand_[0].nstub=int( phcand_[0].h_reg[20]);
      if (phcand_[0].h_reg[60+6]<1.7) continue;
      if ( phcand_[0].h_reg[20]<=0) continue;
      
      if ( phcand_[0].h_reg[70+9]<1.5) continue; //at least 2 Z points
      t.z0=-phcand_[0].h_reg[70+1]/phcand_[0].h_reg[70+0];
      t.eta=phcand_[0].h_reg[70+8];
      if ( fabs(t.z0)>30.) continue;
      
	      
      float theta=GET_THETA_VALUE(ph_,ith);
      float r=GET_R_VALUE(ph_,ir);
      
      double a=-1./tan(theta);
      double b=r/sin(theta);
      
      
      //
      double R=1./2./fabs(r);
      double xi=-a/2./b;
      double yi=1./2./b;
      double g_pt=0.3*3.8*R/100.;
      //g_phi=atan(a);
      double g_phi=theta-PI/2.;
      if (g_phi<0) g_phi+=2*PI;
      CudaHough::Convert(theta,r,&t);
      t.nhits=(pattern>>20)&0x3FF;
      t.theta=theta;
      t.r=r;
      
      t.pt=phcand_[0].h_reg[60+6];
      t.phi=phcand_[0].h_reg[60+2];
      t.nhits=(pattern>>20)&0x3FF;
      getChi2(&phcand_[0],endcap);
      t.chi2=phcand_[0].h_reg[80];
      t.chi2z=phcand_[0].h_reg[81];
      t.layers.clear();
      int layers[2048];
      copyLayers(&phcand_[0],layers);
      for (int ist=0;ist<phcand_[0].nstub;ist++)
	t.layers.push_back(layers[ist]);

      theCandidateVector_.push_back(t);
      
      
      
    }

		 









		  
		
		  
  //  printf("Fin du GPU %ld \n",	theCandidateVector_.size() );



} 
 
void CudaHough::Convert(double theta,double r,mctrack_t *m)
{
  double a=-1./tan(theta);
  double b=r/sin(theta);
		
		
  //
  double R=1./2./fabs(r);
  double xi=-a/2./b;
  double yi=1./2./b;
  //printf(" From r=%f theta=%f a=%f b=%f  R= %f  => Pt=%f GeV/c  Phi0=%f \n",r,theta,a,b,R,0.3*3.8*R/100.,atan(a));
  m->pt=0.3*3.8*R/100.;
  m->phi=atan(a);
  m->phi=theta-PI/2.;
  if (m->phi<0) m->phi+=2*PI;
  m->rho0=sqrt(fabs(R*R-xi*xi-yi*yi));


 
}
void CudaHough::ComputeTracklet(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer)
{int isect=isel;
 DCHistogramHandler* rh=DCHistogramHandler::instance();
  std::stringstream s;
  s<<"/sector"<<isect<<"/";
  TH1* hd3=rh->GetTH1(s.str()+"dist3");
  TH1* hdr3=rh->GetTH1(s.str()+"distr3");
  TH1* hd4=rh->GetTH1(s.str()+"dist4");
  TH1* hdr4=rh->GetTH1(s.str()+"distr4");
  TH1* hphi5=rh->GetTH1("phi5");
  TH1* hc2=rh->GetTH1(s.str()+"chi2");
  TH1* hc2r=rh->GetTH1(s.str()+"chi2r");
  TH1* hr5=rh->GetTH1(s.str()+"chi2");
  TH2* hpthi=rh->GetTH2(s.str()+"pthi");
  if (hd3==0)
    {
      hd3=rh->BookTH1(s.str()+"dist3",2000,-1.,1.);
      hdr3=rh->BookTH1(s.str()+"distr3",200,-5.,5.);
      hd4=rh->BookTH1(s.str()+"dist4",200,-1.,1.);
      hdr4=rh->BookTH1(s.str()+"distr4",200,-10.,10.);
      hphi5=rh->BookTH1("phi5",200,-PI,PI);
      hr5=rh->BookTH1(s.str()+"r5",200,0.,200.);
      hc2=rh->BookTH1(s.str()+"chi2",2000,0.,2.);
      hc2r=rh->BookTH1(s.str()+"chi2r",2000,0.,2.);
      hpthi=rh->BookTH2(s.str()+"pthi",200,0.,20.,200,0.,2*PI);

    }
  


  theNStub_=nstub;
  theX_=x;
  theY_=y;
  theZ_=z;  
  theLayer_=layer;
  memcpy(et_.host_->x_,x,nstub*sizeof(float));
  memcpy(et_.host_->y_,y,nstub*sizeof(float));
  memcpy(et_.host_->z_,z,nstub*sizeof(float));
  memcpy(et_.host_->lay_,layer,nstub*sizeof(uint32_t));
  et_.host_->nstub_=nstub;
  et_.host_->sector_=isel;
  theCandidateVector_.clear();
  // Initialisation depending on sector 
  et_.host_->barrel_=isel>=16 && isel<40;
  et_.host_->inter_=(isel>=8 &&isel<16)||(isel>=40&&isel<48);
  et_.host_->endcap_=(isel<8)||(isel>=48);
  //copyFromHost(&et_);
  et_.host_->ntkl_=0;
  fillDevice(&et_);
  combineLayer(&et_);
  computeTklet(&et_);
  

  addLayer(&et_);
  computeTklet(&et_);
  copyToHost(&et_);
  theCandidateVector_.clear();
  uint32_t ng=0;
  for (int it=0;it<et_.host_->ntkl_;it++)
    {
      ctklet* tk=&(et_.host_->cand_[it]);
      //printf("\t %d %d %x %f %f \n",it,tk->ok_,tk->pattern_,tk->pt_,tk->z0_); 
      if (et_.host_->barrel_ && tk->nxy_<=4) continue;
      if (et_.host_->endcap_ && tk->nxy_<=3) continue;
      if (et_.host_->inter_ && tk->nxy_<=3) continue;
      ctklet tkext;
      memcpy(&tkext,tk,sizeof(ctklet));

      //hpthi->Fill(fabs(tk->pt_),tk->phi_);
      if (fabs(tk->z0_)>20.) continue;
      if (fabs(tk->pt_)<1.8) continue;
      hc2->Fill(TMath::Prob(tk->chi2_,(tk->nxy_-2)));
      if ( TMath::Prob(tk->chi2_,tk->nxy_-2)<1E-3) continue;

      hc2r->Fill(TMath::Prob(tk->chi2r_,tk->nzr_-2));
      if (tk->nzr_>2 && (et_.host_->endcap_ ) &&  TMath::Prob(tk->chi2r_,tk->nzr_-2)<5E-3) continue;
      if (et_.host_->inter_ && tk->nzr_==2) continue;
      //regressiontklet(&tkext,&evt);
      ng++;
      mctrack_t t;
      t.z0=tk->z0_;
      t.eta=tk->eta_;

			  
	      
      t.nhits=tk->nhit_;
      t.theta=tk->theta_;
      t.r=tk->R_;

      t.pt=fabs(tk->pt_);
      t.phi=tk->phi_;

      t.layers.clear();
      theCandidateVector_.push_back(t);
    }
  //printf(" Candidat %d  %d\n",et_.host_->ntkl_,ng);
  
}
