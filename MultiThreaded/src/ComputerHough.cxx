#include <stdint.h>
#include <stdio.h> 
#include "ComputerHough.h"
#include "DCHistogramHandler.h"
#include <math.h>
#include <string.h>
#include <iostream>

ComputerHough::ComputerHough(HoughCut* cuts) :theCuts_(cuts)
{
  theNStub_=0;
  theX_=NULL;
  theY_=NULL;
  theZ_=NULL;
  theLayer_=NULL;

  createHoughCPU(&ph_,768,3072,768);
  for (int i=0;i<96;i++)
    createHoughCPU(&phcand_[i]);
  for (int i=0;i<64;i++)
    createHoughCPU(&phrcand_[i]);


}
void ComputerHough::DefaultCuts()
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
void ComputerHough::Compute(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer)
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
  
  initialiseHoughCPU(&ph_,theNStub_,ntheta,nrho,thmin,thmax,rhmin,rhmax);
  // Rough process
  fillConformalHoughCPU(&ph_,theX_,theY_,theZ_);
  fillLayerHoughCPU(&ph_,theLayer_);
		  //clearHough(&ph);
  processHoughCPU(&ph_,theCuts_->NStubLow,theCuts_->NLayerRow,0,endcap);
  //printf("SECTOR %d gives %d candidates Max val %d STubs %d\n",isel,ph_.h_cand[0],ph_.max_val,ph_.nstub);
  // Precise HT filling
  uint32_t nc=(int)ph_.h_cand[0];
  if (nc>96) nc=96;
  for (int ic=0;ic<nc;ic++)
    {
      clearHoughCPU(&phcand_[ic]);
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
      
      initialiseHoughCPU(&phcand_[ic],theNStub_,nbinf,nbinr,tmi,tma,rmi,rma);	    

      copyPositionHoughCPU(&ph_,pattern,&phcand_[ic],0,false,endcap);
    }
		  

		
  //Precise HT processing
		 
  for (int ic=0;ic<nc;ic++)
    {
      if (phcand_[ic].h_reg[20]>0)
	{
	  phcand_[ic].nstub=int( phcand_[ic].h_reg[20]);
	  processHoughCPU(&phcand_[ic],theCuts_->NStubHigh,theCuts_->NLayerHigh,0,endcap);
	  
	}
    }

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
	      initialiseHoughCPU(&phrcand_[ici],theNStub_,32,32,-PI/2,PI/2,-150.,150.);
	      copyPositionHoughCPU(&phcand_[ic],patterni,&phrcand_[ici],1,true,endcap);
	      phrcand_[ici].nstub=int( phrcand_[ici].h_reg[20]);
	      if (phrcand_[ici].h_reg[60+6]<1.7) continue;
	      if ( phrcand_[ici].h_reg[20]<=0) continue;
			      
	      if ( phrcand_[ici].h_reg[70+9]<1.5) continue; //at least 2 Z points
	      //@@@@@ un essai
	      t.z0=-phrcand_[ici].h_reg[70+1]/phrcand_[ici].h_reg[70+0];
	      t.eta=phrcand_[ici].h_reg[70+8];
	      if ( fabs(t.z0)>30.) continue;
			  
	      
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
	      ComputerHough::Convert(theta,r,&t);
	      t.nhits=(patterni>>20)&0x3FF;
	      t.theta=theta;
	      t.r=r;

	      t.pt=phrcand_[ici].h_reg[60+6];
	      t.phi=phrcand_[ici].h_reg[60+2];
	      t.nhits=(patterni>>20)&0x3FF;
	      t.layers.clear();
	      for (int ist=0;ist<phrcand_[ici].nstub;ist++)
		t.layers.push_back(phrcand_[ici].d_layer[ist]);
	      theCandidateVector_.push_back(t);

			      
	    }
	}
    }
		 









		  
		
		  
  //  printf("Fin du CPU %ld \n",	theCandidateVector_.size() );



}
void ComputerHough::ComputeOneShot(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer)
{
  //ComputeTracklet(isel,nstub,x,y,z,layer);
  //return;
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

  initialiseHoughCPU(&ph_,theNStub_,ntheta,nrho,thmin,thmax,rhmin,rhmax);
  // Rough process
  fillConformalHoughCPU(&ph_,theX_,theY_,theZ_);
  fillLayerHoughCPU(&ph_,theLayer_);
		  //clearHough(&ph);
  processHoughCPU(&ph_,theCuts_->NStubLow,theCuts_->NLayerRow,0,endcap);
  //printf("SECTOR %d gives %d candidates Max val %d STubs %d\n",isel,ph_.h_cand[0],ph_.max_val,ph_.nstub);
  // Precise HT filling
  uint32_t nc=(int)ph_.h_cand[0];
  if (nc>512) nc=512;
  clearHoughCPU(&phcand_[0]);
  

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
      initialiseHoughCPU(&phcand_[0],theNStub_,32,32,-PI/2,PI/2,-150.,150.);
      copyPositionHoughCPU(&ph_,pattern,&phcand_[0],1,true,endcap);
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
	//printf("PTTTT %f %f \n",g_pt,5.7E-3*sqrt(a*a+1)/b);
      double g_phi=theta-PI/2.;
      if (g_phi<0) g_phi+=2*PI;
      ComputerHough::Convert(theta,r,&t);
      t.nhits=(pattern>>20)&0x3FF;
      t.theta=theta;
      t.r=r;
      
      t.pt=phcand_[0].h_reg[60+6];
      t.phi=phcand_[0].h_reg[60+2];
      t.nhits=(pattern>>20)&0x3FF;
      t.layers.clear();
      for (int ist=0;ist<phcand_[0].nstub;ist++)
	t.layers.push_back(phcand_[0].d_layer[ist]);
      theCandidateVector_.push_back(t);
      
      
      
    }
}
 
void ComputerHough::Convert(double theta,double r,mctrack_t *m)
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
#include <sstream>
#include <iostream>
typedef struct {

  int32_t pattern_;
  int32_t nhit_;
  uint16_t idx_[32];
  double sumx_,sumx2_,sumxy_,sumy_,nxy_;
  double sumz_,sumz2_,sumzr_,sumr_,nzr_;
  double ax_,bx_,phi_,theta_,R_,pt_;
  double ar_,br_,eta_,z0_,chi2_,chi2r_;

} tklet;

#define GPU_MAX_TKLET 8192
typedef struct
{
  int sector_,nstub_;
  float x_[GPU_MAX_STUB],y_[GPU_MAX_STUB];
  float z_[GPU_MAX_STUB],r_[GPU_MAX_STUB];
  float xp_[GPU_MAX_STUB],yp_[GPU_MAX_STUB];
  uint32_t lay_[GPU_MAX_STUB];
  int ntkl_;
  tklet cand_[GPU_MAX_TKLET];
  bool barrel_,inter_,endcap_;
} tklevent;

void combineLayer(tklevent* e,uint32_t l1,uint32_t l2)
{
  for (uint32_t i=0;i<e->nstub_;i++)
    {
      if ((e->lay_[i]&0xFFFF)!=l1) continue;
      if (((e->lay_[i]>>16)&0x3)==0) continue;
      for (uint32_t j=0;j<e->nstub_;j++)
	{
	  if ((e->lay_[j]&0xFFFF)!=l2) continue;
	  if (((e->lay_[j]>>16)&0x3)==0) continue;
	   // calcul de a b en rphi
	  double a =(e->yp_[j]-e->yp_[i])/(e->xp_[j]-e->xp_[i]);
	  //if ((isel%4==0) && (a<-0.25 || a>1.55)) continue;
	  double b =e->yp_[j]-a*e->xp_[j];
	  double pt=5.7E-3*sqrt((a*a+1)/b/b);
	  if (fabs(pt)<1.8) continue;
	  double ar=(e->r_[j]-e->r_[i])/(e->z_[j]-e->z_[i]);
	  double br=e->r_[j]-ar*e->z_[j];
	  double zi=-br/ar;
	  if (fabs(zi)>20) continue;
	  tklet* tkl=&e->cand_[e->ntkl_];
	  //printf("%x %x %d %d \n",e->lay_[i],e->lay_[j],(e->lay_[i]>>16)&0x3,(e->lay_[j]>>16)&0x3);
	  tkl->nhit_=0; 
	  tkl->pattern_=0; 
	  tkl->pattern_ |=(1<<l1);
	  tkl->pattern_ |=(1<<l2);
	  tkl->idx_[tkl->nhit_++]=i;
	  tkl->idx_[tkl->nhit_++]=j;
	  e->ntkl_++;
	}
    }
} 

void  regressiontklet(tklet* t,tklevent* e)
{
  DCHistogramHandler* rh=DCHistogramHandler::instance();
  for (int ilay=5;ilay<24;ilay++)
    {
      if ((t->pattern_ & (1<<ilay))==0) continue;
      memset(&t->sumx_,0,20*sizeof(double));
  //printf("clearing %d %d \n",it,t->nhit_);
      int idx;
      for (int ih=0;ih<t->nhit_;ih++)
	{
	  idx=t->idx_[ih];
	  int layerh=(e->lay_[idx]&0xFFFF);
	  if (layerh==ilay) continue;
	  t->sumx_+=e->xp_[idx];
	  t->sumx2_+=e->xp_[idx]*e->xp_[idx];
	  t->sumxy_+=e->xp_[idx]*e->yp_[idx];
	  t->sumy_+=e->yp_[idx];
	  t->nxy_+=1.;
	  if (((e->lay_[idx]>>16)&0x3)==0) continue;
	  t->sumz_+=e->z_[idx];
	  t->sumz2_+=e->z_[idx]*e->z_[idx];
	  t->sumzr_+=e->z_[idx]*e->r_[idx];
	  t->sumr_+=e->r_[idx];
	  t->nzr_+=1.;
      
	}
      if (t->nxy_<2) continue; 
      double s2z = t->sumx2_/t->nxy_-(t->sumx_/t->nxy_)*(t->sumx_/t->nxy_);
      double szx = t->sumxy_/t->nxy_-(t->sumx_/t->nxy_)*(t->sumy_/t->nxy_);
  
      t->ax_ = szx/s2z;
      t->bx_=(t->sumy_/t->nxy_)-t->ax_*(t->sumx_/t->nxy_);
      // Now extrapolate to layer ilay
      //printf(" extrapolating %d \n",ilay);
      std::stringstream s;
      s<<"/sector"<<e->sector_<<"/layer"<<ilay;

      TH1* hdistx=rh->GetTH1(s.str()+"/distxy");
      TH1* hdxymin=rh->GetTH1(s.str()+"/dxymin");
      if (hdistx==0)
	{
	  hdistx=rh->BookTH1(s.str()+"/distxy",300,-0.5E-3,0.5E-3);
	  hdxymin=rh->BookTH1(s.str()+"/dxymin",300,-0.5E-3,0.5E-3);
	  printf("Booking %s histo\n",s.str().c_str());
	}
      double dmin=9E10,admin=9E10;
      for (int is=0;is<e->nstub_;is++)
	{
	  //int is=t->idx_[ih];

	  int layerh=(e->lay_[is]&0xFFFF);
	  if (layerh!=ilay) continue;
	  double r2=e->r_[is]*e->r_[is];
	  r2=1.;
	  double distx=r2*(t->ax_*e->xp_[is]+t->bx_-e->yp_[is])/sqrt(1+t->ax_*t->ax_);
	  //	  if (fabs(distx)<1E-3)
	  //printf("is %d %f r2 %f %f %f %x pattern %x  \n",is,e->r_[is],distx,r2,t->ax_,e->lay_[is],t->pattern_);
	  

	  if (fabs(distx)<admin)
	    {
	      admin=fabs(distx);
	      dmin=distx;
		
	    }
	  hdistx->Fill(distx);
	}
      hdxymin->Fill(dmin);
    }
}

uint32_t computeTklets(tklevent* e)
{
  uint32_t ngood=0;
  //printf ("Number of tklet %d\n",e->ntkl_);
  for (int it=0;it<e->ntkl_;it++)
    {
      tklet* t=&(e->cand_[it]);
      if (t->nhit_<=0) continue;
      t->pattern_=0;
      
      memset(&t->sumx_,0,20*sizeof(double));
      //printf("clearing %d %d \n",it,t->nhit_);
      int idx;
      for (int ih=0;ih<t->nhit_;ih++)
	{
	  idx=t->idx_[ih];
	  int layerh=(e->lay_[idx]&0xFFFF);
	  t->pattern_ |=(1<<layerh);
	  t->sumx_+=e->xp_[idx];
	  t->sumx2_+=e->xp_[idx]*e->xp_[idx];
	  t->sumxy_+=e->xp_[idx]*e->yp_[idx];
	  t->sumy_+=e->yp_[idx];
	  t->nxy_+=1.;
	  if (((e->lay_[idx]>>16)&0x3)==0) continue;
	  t->sumz_+=e->z_[idx];
	  t->sumz2_+=e->z_[idx]*e->z_[idx];
	  t->sumzr_+=e->z_[idx]*e->r_[idx];
	  t->sumr_+=e->r_[idx];
	  t->nzr_+=1.;

	}
      if (t->nzr_<2 || t->nxy_<2) 
	{
	  if (t->nhit_>0) t->nhit_*=-1;
	  continue;
	}
      double s2z = t->sumx2_/t->nxy_-(t->sumx_/t->nxy_)*(t->sumx_/t->nxy_);
      double szx = t->sumxy_/t->nxy_-(t->sumx_/t->nxy_)*(t->sumy_/t->nxy_);
  
      t->ax_ = szx/s2z;
      t->bx_=(t->sumy_/t->nxy_)-t->ax_*(t->sumx_/t->nxy_);

     
      t->phi_=atan(t->ax_);
      if (t->phi_<0) t->phi_+=2*PI;
      float xp1=e->xp_[t->idx_[0]];
      float yp1=e->yp_[t->idx_[0]];
      if (xp1>0 && yp1>0 && t->phi_>PI) t->phi_-=PI;
      if (xp1<0 && yp1>0 && t->phi_>PI) t->phi_-=PI;
      if (xp1<0 && yp1<0 && t->phi_<PI) t->phi_+=PI;
      if (xp1>0 && yp1<0 && t->phi_<PI) t->phi_+=PI;
      t->theta_=atan(-1./t->ax_);
      t->pt_=5.7E-3*sqrt((t->ax_*t->ax_+1)/t->bx_/t->bx_);
      //      printf("%f %f  \n",t->ax_,t->bx_);
      if (fabs(t->pt_)<1.8) {if (t->nhit_>0) t->nhit_*=-1;continue;}
      s2z = t->sumz2_/t->nzr_-(t->sumz_/t->nzr_)*(t->sumz_/t->nzr_);
      szx = t->sumzr_/t->nzr_-(t->sumz_/t->nzr_)*(t->sumr_/t->nzr_);

      t->ar_= szx/s2z;
      t->br_=(t->sumr_/t->nzr_)-t->ar_*(t->sumz_/t->nzr_);
      t->z0_=-t->br_/t->ar_;
      if (fabs(t->z0_)>20) {if (t->nhit_>0) t->nhit_*=-1;continue;}
      t->eta_=-log(fabs(tan(atan( t->ar_)/2)));
      t->eta_=-log(fabs((1+sqrt(1+t->ar_*t->ar_))/t->ar_) );
      if (t->ar_>0) t->eta_=-1*t->eta_;
      if (t->nxy_>4) ngood++;
      
      t->chi2_=0;
      t->chi2r_=0;
      for (int ih=0;ih<t->nhit_;ih++)
	{
	  int is=t->idx_[ih];
	  int idx=is;
	  //double delta=(t->ax_*e->xp_[idx]+t->bx_-e->yp_[idx])*e->r_[idx]*e->r_[idx]/0.04;
	  double r2=e->r_[is]*e->r_[is];
	  double distx=r2*(t->ax_*e->xp_[is]+t->bx_-e->yp_[is])/sqrt(1+t->ax_*t->ax_)/0.03;
	  distx=(t->ax_*e->xp_[is]+t->bx_-e->yp_[is])/sqrt(1+t->ax_*t->ax_)/8.4E-6;
	  t->chi2_+=distx*distx;
	  if (((e->lay_[idx]>>16)&0x3)==0) continue;
	  double deltar=(t->ar_*e->z_[idx]+t->br_-e->r_[idx])/sqrt(1+t->ar_*t->ar_)/0.06;
	  t->chi2r_+=deltar*deltar;

	}
    
    }
  return ngood;
}
void addLayer(tklevent* e,uint32_t l1,TH1* hdist, TH1* hdistr)
{
   uint32_t ngood=0;
  for (int it=0;it<e->ntkl_;it++)
    {
      tklet* t=&(e->cand_[it]);
      if (t->nhit_<=0) continue;
      if ((t->pattern_&(1<<l1))!=0) continue; // layer already included

      for (int is=0;is<e->nstub_;is++)
	{
	  if ((e->lay_[is]&0xFFFF)!=l1) continue;
	  double r2=e->r_[is]*e->r_[is];
	  r2=1.;
	  double distx=r2*(t->ax_*e->xp_[is]+t->bx_-e->yp_[is])/sqrt(1+t->ax_*t->ax_);
	  //	  if (fabs(distx)<1E-3)
	  //printf("is %d %f r2 %f %f %f %x pattern %x  \n",is,e->r_[is],distx,r2,t->ax_,e->lay_[is],t->pattern_);

	  if (hdist!=0) hdist->Fill(distx*e->r_[is]*e->r_[is]);
	  //if (fabs(distx)>0.3) continue;
	  double cut=0.15;
	  cut=3.5E-5;
	  if (l1<=10) cut /=2;
	  if (e->barrel_ && fabs(distx)>cut) continue;
	  if (e->inter_ && fabs(distx)>cut) continue;
	  if (e->endcap_ && fabs(distx)>cut) continue;
	  if (((e->lay_[is]>>16)&0x3)!=0)
	    {
	      double distr=(t->ar_*e->z_[is]+t->br_-e->r_[is])/sqrt(1+t->ar_*t->ar_);
	      if (hdistr!=0) hdistr->Fill(distr);
	      if (fabs(distr)>0.6) continue;
	      //  if (e->barrel_ && fabs(distr)>0.45) continue;
	      // if (e->inter_ && fabs(distr)>0.2) continue;
	      //if (e->endcap_ && fabs(distr)>0.2) continue;
	    }
	  t->idx_[t->nhit_++]=is;
	  t->pattern_|=(1<<l1);


	  break;
	}
    }
}
static tklevent evt;
void ComputerHough::ComputeTracklet(uint32_t isel,uint32_t nstub,float* x,float* y,float* z,uint32_t* layer,int32_t* flayer)
{
 int isect=isel;
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
  

  theCandidateVector_.clear();
 
 theNStub_=nstub;
  theX_=x;
  theY_=y;
  theZ_=z;  
  theLayer_=layer;
 
  bool barrel=isect>=16 && isect<40;
  bool inter=(isect>=8 &&isect<16)||(isect>=40&&isect<48);
  bool endcap=(isect<8)||(isect>=48);
  evt.barrel_=barrel;
  evt.inter_=inter;
  evt.endcap_=endcap;
  float thmin=-PI/2,thmax=PI/2;
  float rhmin=theCuts_->RhoMin,rhmax=theCuts_->RhoMax;
  //printf("On appelle le GPU %d \n",theNStub_);
  int ntheta=160;
  int nrho=theCuts_->NRho;//8//12//192;
  //initialiseHough(&ph,gpu_nstub,ntheta,nrho,-PI/2,PI/2,-0.06,0.06);
  if (barrel || inter || endcap)
    {
      ntheta=theCuts_->NTheta;//64;
      if (isect%4==0) thmin=1.32;
      if (isect%4==1) thmin=-1.04;
      if (isect%4==2) thmin=-0.24;
      if (isect%4==3) thmin=0.51;
      thmax=thmin+1.25;
    }
  evt.sector_=isect;
  evt.nstub_ = nstub;
  memcpy(evt.x_,x,nstub*sizeof(float));
  memcpy(evt.y_,y,nstub*sizeof(float));
  memcpy(evt.z_,z,nstub*sizeof(float));
  memcpy(evt.lay_,layer,nstub*sizeof(uint32_t));
  memset(&evt.cand_,0,GPU_MAX_TKLET*sizeof(tklet));
  for (uint16_t is=0;is<nstub;is++)
    {
      double r2=(x[is]*x[is]+y[is]*y[is]);
      evt.xp_[is]=x[is]/r2;
      evt.yp_[is]=y[is]/r2;
      
      evt.r_[is]=sqrt(r2);
      if (evt.r_[is]<1.)
	printf("%d is %f %f %f %f \n",is,r2,evt.r_[is],x[is],y[is]);
    }
  evt.ntkl_=0;
  if (inter || barrel)
    {
      combineLayer(&evt,5,6);
      //printf("5-6 %d \n",evt.ntkl_);
      combineLayer(&evt,5,7);
      //printf("5-7 %d \n",evt.ntkl_);
      combineLayer(&evt,6,7);
    }
  else
    if (isect<8)
      {
	combineLayer(&evt,5,18);
      //printf("5-18 %d \n",evt.ntkl_);
	combineLayer(&evt,5,19);
      //printf("5-19 %d \n",evt.ntkl_);
	combineLayer(&evt,18,19);

      }
    else
      if (isect>=48)
	{
	  combineLayer(&evt,5,11);
	  //printf("5-11 %d \n",evt.ntkl_);
	  combineLayer(&evt,5,12);
	  // printf("5-12 %d \n",evt.ntkl_);
	  combineLayer(&evt,11,12);

	}
  //printf("Number of tracklet %d \n",evt.ntkl_);
  uint32_t ng= computeTklets(&evt);
  //printf("nb 2 poins %d \n",ng);
  for (int il=24;il>=5;il--)
    { 
      addLayer(&evt,il,hd3,hdr3);
      //ng= computeTklets(&evt);
      
    }
  ng= computeTklets(&evt);
  ng=0;
  theCandidateVector_.clear();
  for (int it=0;it<evt.ntkl_;it++)
    {
      tklet* tk=&(evt.cand_[it]);
      if (barrel && tk->nxy_<=4) continue;
      if (endcap && tk->nxy_<=3) continue;
      if (inter && tk->nxy_<=3) continue;
      tklet tkext;
      memcpy(&tkext,tk,sizeof(tklet));

      hpthi->Fill(fabs(tk->pt_),tk->phi_);
      if (fabs(tk->z0_)>20.) continue;
      if (fabs(tk->pt_)<1.8) continue;
      hc2->Fill(TMath::Prob(tk->chi2_,(tk->nxy_-2)));
      if ( TMath::Prob(tk->chi2_,tk->nxy_-2)<1E-3) continue;

      hc2r->Fill(TMath::Prob(tk->chi2r_,tk->nzr_-2));
      if (tk->nzr_>2 && (endcap ) &&  TMath::Prob(tk->chi2r_,tk->nzr_-2)<5E-3) continue;
      if (inter && tk->nzr_==2) continue;
      //regressiontklet(&tkext,&evt);
      ng++;
      mctrack_t t;
      t.z0=tk->z0_;
      t.eta=tk->eta_;
      t.matches=tk->pattern_;
			  
	      
      t.nhits=tk->nhit_;
      t.theta=tk->theta_;
      t.r=tk->R_;

      t.pt=fabs(tk->pt_);
      t.phi=tk->phi_;

      t.layers.clear();
      theCandidateVector_.push_back(t);
    }
  //printf("Candidats %d \n",ng);
  return;
  
}
 
