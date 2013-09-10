#include <HoughCartesian.h>
#include <TLine.h>
#include <TEllipse.h>

#include <sstream>
HoughCartesian::HoughCartesian(double thmin,double thmax,double rmin,double rmax,uint32_t nbintheta,uint32_t nbinr) : theThetaMin_(thmin),theThetaMax_(thmax),theRMin_(rmin),theRMax_(rmax),theNbinTheta_(nbintheta),theNbinR_(nbinr)
{
	this->clear();
	theThetaBin_=(theThetaMax_-theThetaMin_)*1./theNbinTheta_;
	theRBin_=(theRMax_-theRMin_)*1./theNbinR_;
	

}
HoughCartesian::~HoughCartesian(){;}
uint32_t HoughCartesian::getVoteMax(){return theVoteMax_;}
void HoughCartesian::fill(double z, double x)
{
	
	//theVoteMax_=0;
	theX_.push_back(z);
	theY_.push_back(x);
	double zp=1*z/(z*z+x*x);
	double xp=1*x/(z*z+x*x);
	
	//printf("x %f y %f \n",zp,xp);
	for (uint32_t i=0;i<theNbinTheta_;i++)
				{
				  double r = xp-getTheta(i)*zp;
					//printf("R %f Cos %f \n",r,theCos_[i]);
					if (r>theRMax_ || r<theRMin_) continue;
					uint32_t ir=int(floor((r-theRMin_)/theRBin_));
					theHoughImage_[i][ir]+=1;
					if (theHoughImage_[i][ir]>theVoteMax_) theVoteMax_=theHoughImage_[i][ir];
				}
		
}
void HoughCartesian::clear() 
{
	memset(theHoughImage_,0,1024*1024*sizeof(uint16_t));
	theX_.clear();
	theY_.clear();
	theVoteMax_=0;
}
double HoughCartesian::getTheta(int32_t i) {return theThetaMin_+i*theThetaBin_;}
double HoughCartesian::getR(int32_t i) {return theRMin_+i*theRBin_;}
uint16_t HoughCartesian::getValue(uint32_t i,uint32_t j){return theHoughImage_[i][j];}
void HoughCartesian::findMaxima(std::vector< std::pair<uint32_t,uint32_t> >& maxval,uint32_t cut)
{
		printf("The Vote maximum is %d \n",theVoteMax_);
#ifdef VOTEBIN
	if (theVoteMax_ < 3) return;
	double theVoteBin_=theVoteMax_*1./(100);
	printf("the Vote bin %f \n",theVoteBin_);
	uint16_t theVotes_[100];
	memset(theVotes_,0,100*sizeof(uint16_t));
	for (uint32_t ith=0;ith<theNbinTheta_;ith++)
		for (uint32_t ir=0;ir<theNbinR_;ir++)
			for (uint32_t iv=0;iv<100;iv++) if (theHoughImage_[ith][ir]>iv*theVoteBin_) theVotes_[iv]++;
	
	uint32_t ivmin=0;
	for (uint32_t iv=100-1;iv>=0;iv--) if (theVotes_[iv]>cut) {ivmin=iv+1;break;}
	if (ivmin == 100) ivmin--;
	double vcut=ivmin*theVoteBin_;
#else
	double vcut=theVoteMax_-1;
	if (vcut>12) vcut=12;
#endif
	printf("The vote cut is %d \n",cut);
	for (uint32_t ith=0;ith<theNbinTheta_;ith++)
		for (uint32_t ir=0;ir<theNbinR_;ir++)
			{
				uint32_t count=theHoughImage_[ith][ir];
				if (count<cut) continue;
				bool notmax=false;
					for (int ithb=-1;ithb<=1;ithb++)
						for (int irb=-1;irb<=1;irb++)
						{
							if (ithb==0 && irb==0) continue;
							if ((ith+ithb)<0 || (ith+ithb)>theNbinTheta_) continue;
							if ((ir+irb)<0 || (ir+irb)>theNbinR_) continue;
							notmax= notmax || (count<theHoughImage_[ith+ithb][ir+irb]);
							}
					if (notmax) continue;
					for (std::vector< std::pair<uint32_t,uint32_t> >::iterator it=maxval.begin();it!=maxval.end();it++)
						if (TMath::Abs((int)(ith-it->first))<2 &&TMath::Abs((int)(ir-it->second))<2)
							{ notmax=true;break;}
					if (notmax) continue;
					std::pair <uint32_t,uint32_t> p(ith,ir);
					maxval.push_back(p);
				}
	
}
void HoughCartesian::findMaximumBins(std::vector< std::pair<double,double> >& maxval,uint32_t cut)
{



int32_t** matrix= new int*[theNbinTheta_];
for(int i = 0; i < theNbinTheta_; ++i) {
    matrix[i] = new int[theNbinR_];
}

 int32_t** matlap= new int*[theNbinTheta_];
for(int i = 0; i < theNbinTheta_; ++i) {
    matlap[i] = new int[theNbinR_];
} 



for (uint32_t ith=0;ith<theNbinTheta_;ith++)
		for (uint32_t ir=0;ir<theNbinR_;ir++)
			{matrix[ith][ir]=0;matlap[ith][ir]=0;}
			

for (uint32_t ith=0;ith<theNbinTheta_;ith++)
		for (uint32_t ir=0;ir<theNbinR_;ir++)
			{
				uint32_t count=theHoughImage_[ith][ir];
				if (count<cut) continue;
				bool notmax=false;
					for (int ithb=-1;ithb<=1;ithb++)
						for (int irb=-1;irb<=1;irb++)
						{
							if (ithb==0 && irb==0) continue;
							if ((ith+ithb)<0 || (ith+ithb)>theNbinTheta_) continue;
							if ((ir+irb)<0 || (ir+irb)>theNbinR_) continue;
							notmax= notmax || (count<theHoughImage_[ith+ithb][ir+irb]);
							}
					if (notmax) continue;
					matrix[ith][ir]=1;
					
				}
	
	// Now merge adjacent bins
		int weight[3][3]={{1,1,1},{1,-8,1},{1,1,1}};
		for (uint32_t ith=1;ith<theNbinTheta_-1;ith++)
		for (uint32_t ir=1;ir<theNbinR_-1;ir++)
			{
			for (int j = -1; j <= 1; j++) {
				for (int i = -1; i <= 1; i++) {
					matlap[ith][ir] += weight[j + 1][i + 1] * matrix[ith + j][ir + i];
				}
				}
				if (matlap[ith][ir]==-8 && matrix[ith][ir]!=0) //isolated maximum
					{
					  //std::cout<<"standalone " <<ith<<":"<<ir<<std::endl;
						std::pair<double,double> pf((ith+0.5)*theThetaBin_+theThetaMin_,(ir+0.5)*theRBin_+theRMin_);
						maxval.push_back(pf);
					}
			}
	// Make the Mean of adjacent bins
	for (uint32_t ith=1;ith<theNbinTheta_-1;ith++)
		for (uint32_t ir=1;ir<theNbinR_-1;ir++)
			{
				if (matrix[ith][ir]==0 || matlap[ith][ir] == -8) continue;
				//std::cout<<" central= "<<ith<<":"<<ir<<"="<<matlap[ith][ir]<<std::endl;
				double theta=0,r=0,n=0;
				for (int j = -1; j <= 1; j++) {
				for (int i = -1; i <= 1; i++) {
					if (matlap[ith+j][ir+i] > -8 && matrix[ith+j][ir+i]!=0  )
						{
							//std::cout<<" adj= "<<ith+j<<":"<<ir+i<<"="<<matlap[ith+j][ir+i]<<std::endl;

							theta+=(ith+j)+0.5;
							r+=(ir+i)+0.5;
							n+=1;
							matrix[ith+j][ir+i]=0;
						}
				}
				}
				theta/=n;
				r/=n;
				//std::cout<<"Nearby " <<theta<<":"<<r<<std::endl;
				std::pair<double,double> pf(theta*theThetaBin_+theThetaMin_,r*theRBin_+theRMin_);
				maxval.push_back(pf);
			}

for(int i = 0; i < theNbinTheta_; ++i) 
	{
    delete [] matrix[i];
  }
	delete [] matrix;
	for(int i = 0; i < theNbinTheta_; ++i) {
    delete [] matlap[i];
}
	delete [] matlap; 
}


void HoughCartesian::findMaxima(std::vector< std::pair<double,double> >& maxval,uint32_t cut)
{
		
	for (uint32_t ith=0;ith<theNbinTheta_;ith++)
	for (uint32_t ir=0;ir<theNbinR_;ir++)
	if (theHoughImage_[ith][ir]>=cut)
	{
		bool notmax=false;
					for (int ithb=-1;ithb<=1;ithb++)
						for (int irb=-1;irb<=1;irb++)
						{
							if (ithb==0 && irb==0) continue;
							if ((ith+ithb)<0 || (ith+ithb)>theNbinTheta_) continue;
							if ((ir+irb)<0 || (ir+irb)>theNbinR_) continue;
							notmax= notmax || (theHoughImage_[ith][ir]<theHoughImage_[ith+ithb][ir+irb]);
							}
					if (notmax) continue;
					
		double theta=0,r=0,w=0,nb=0;
		
		for (int ithb=-1;ithb<=1;ithb++)
		for (int irb=-1;irb<=1;irb++)
		{
			if ((ith+ithb)<0 || (ith+ithb)>theNbinTheta_) continue;
			if ((ir+irb)<0 || (ir+irb)>theNbinR_) continue;
			double thetab = (this->getTheta(ith+ithb)+this->getTheta(ith+ithb+1))/2;
			double rb = (this->getR(ir+irb)+this->getR(ir+irb+1))/2;
			double wb=theHoughImage_[ith+ithb][ir+irb]*1.;
			if (ithb!=0 || irb!=0) wb=wb/2.;
			//printf("%d %d %f \n",ithb,irb,wb);
			w+=wb;
			theta+=thetab*wb;
			r+=rb*wb;
			nb+=1;
		}
		printf("\t Candidat %f %f %f \n",w,theta/w,r/w);
		std::pair <double,double> p(theta/w,r/w/1.);
		//std::pair <double,double> p((ith+0.5)*theThetaBin_+theThetaMin_,(ir+0.5)*theRBin_+theRMin_);
		maxval.push_back(p);
	}
	
}

static TCanvas* CanvasHough=NULL;
void HoughCartesian::draw(DCHistogramHandler* h,std::vector< std::pair<uint32_t,uint32_t> > *maxval)
{
if (CanvasHough==NULL)
	{
		CanvasHough=new TCanvas("CanvasHough","hough",800,900);
		CanvasHough->Modified();
		CanvasHough->Draw();
		CanvasHough->Divide(1,2);
		TVirtualPad* pd=CanvasHough->cd(2);
		pd->Divide(2,1);
	}
	CanvasHough->cd();
	std::stringstream htname;
	htname<<"HoughCartesian"<<theNbinTheta_<<"_"<<theNbinR_;
	TH2F* hhtx = (TH2F*) h->GetTH2(htname.str());
	TH2F* hx = (TH2F*) h->GetTH2("LocalInverse");
	TH2F* hox = (TH2F*) h->GetTH2("LocalImage");

	if (hhtx==NULL)
	{
		hhtx =(TH2F*)h->BookTH2(htname.str(),theNbinTheta_,theThetaMin_,theThetaMax_,theNbinR_,theRMin_,theRMax_);
	}
	else
		hhtx->Reset();

	if (hx==NULL)
	{
		hx =(TH2F*)h->BookTH2("LocalInverse",100,0.,0.05,100,0.,0.05);
		hox =(TH2F*)h->BookTH2("LocalImage",800,-50.,150.,800,-50.,150);
		hx->SetMarkerColor(4);
	hx->SetMarkerSize(.2);
	hx->SetMarkerStyle(25);
	hox->SetMarkerColor(4);
	hox->SetMarkerSize(.2);
	hox->SetMarkerStyle(25);
	}
	else
	{
		hx->Reset();
		hox->Reset();
	}
for (int ip=0;ip<theX_.size();ip++)
		{
			double xp=theX_[ip]/(theX_[ip]*theX_[ip]+theY_[ip]*theY_[ip]);
			double yp=theY_[ip]/(theX_[ip]*theX_[ip]+theY_[ip]*theY_[ip]);
			hx->Fill(xp,yp);
			hox->Fill(theX_[ip],theY_[ip]);
		
		}
		
	
	for (uint32_t i=0;i<theNbinTheta_;i++)
	for (uint32_t j=0;j<theNbinR_;j++)
	{
		hhtx->SetBinContent(i+1,j+1,theHoughImage_[i][j]*1.);

		
	}
	
	CanvasHough->cd(1);
	hhtx->Draw("COLZ");
	//hw->Draw();
	TVirtualPad* pd=CanvasHough->cd(2);
	pd->cd(1);
	

	hx->Draw("p");
	std::vector<TLine*> vline;vline.clear();
	if (maxval!=NULL)
	{
		for (std::vector < std::pair<uint32_t,uint32_t> >::iterator ihb=maxval->begin();ihb<maxval->end();ihb++)
		{
		uint32_t ith=(*ihb).first;
		uint32_t ir=(*ihb).second;
		double theta=0,r=0,w=0,nb=0;
		
		for (int ithb=-1;ithb<=1;ithb++)
		for (int irb=-1;irb<=1;irb++)
		{
		 if ((ith+ithb)<0 || (ith+ithb)>theNbinTheta_) continue;
		 if ((ir+irb)<0 || (ir+irb)>theNbinR_) continue;
		double thetab = (this->getTheta(ith+ithb)+this->getTheta(ith+ithb+1))/2;
		double rb = (this->getR(ir+irb)+this->getR(ir+irb+1))/2;
		double wb=theHoughImage_[ith+ithb][ir+irb]*1.;
		//printf("%d %d %f \n",ithb,irb,wb);
		w+=wb;
		theta+=thetab*wb;
		r+=rb*wb;
		nb+=1;
		}
		r=r/w;
		theta=theta/w;
		double a=-1./tan(theta);
		double b=r/sin(theta);
		
		
		//
		double R=1./2./TMath::Abs(r);
		double xi=-a/2./b;
		double yi=1./2./b;
		printf("%f %f %f %f %f %f Rayon R= %f => %f GeV/c\n",nb,w,r,theta,a,b,R,0.3*3.8*R/100.);
		TLine* l =new TLine(0.,b,0.05,a*0.05+b);
			l->SetLineColor(2);
			l->Draw("SAME");
				vline.push_back(l);
			if (vline.size()>5) break;
		}
	}
	pd->cd(2);
	hox->Draw("p");
	std::vector<TEllipse*> vel;vel.clear();
	if (maxval!=NULL)
	{
		for (std::vector < std::pair<uint32_t,uint32_t> >::iterator ihb=maxval->begin();ihb<maxval->end();ihb++)
		{
		uint32_t ith=(*ihb).first;
		uint32_t ir=(*ihb).second;
		double theta=0,r=0,w=0,nb=0;
		
		for (int ithb=-1;ithb<=1;ithb++)
		for (int irb=-1;irb<=1;irb++)
		{
		 if ((ith+ithb)<0 || (ith+ithb)>theNbinTheta_) continue;
		 if ((ir+irb)<0 || (ir+irb)>theNbinR_) continue;
		double thetab = (this->getTheta(ith+ithb)+this->getTheta(ith+ithb+1))/2;
		double rb = (this->getR(ir+irb)+this->getR(ir+irb+1))/2;
		double wb=theHoughImage_[ith+ithb][ir+irb]*1.;
		//printf("%d %d %f \n",ithb,irb,wb);
		w+=wb;
		theta+=thetab*wb;
		r+=rb*wb;
		nb+=1;
		}
		r=r/w;
		theta=theta/w;
		double a=-1./tan(theta);
		double b=r/sin(theta);
		
		
		//
		double R=1./2./TMath::Abs(r);
		double xi=-a/2./b;
		double yi=1./2./b;
		printf("%f %f %f %f %f %f Rayon R= %f => %f GeV/c\n",nb,w,r,theta,a,b,R,0.3*3.8*R/100.);
		TEllipse* l =new TEllipse(xi,yi,R,R);
			l->SetLineColor(3);
			l->SetFillStyle(0);
			l->Draw("SAME");
				vel.push_back(l);
			if (vel.size()>5) break;
		}
	}
	CanvasHough->Modified();
	CanvasHough->Draw();
	//CanvasHough->WaitPrimitive();

	CanvasHough->Update();

	char c;c=getchar();putchar(c); if (c=='.') exit(0);
	for (std::vector<TLine*>::iterator il=vline.begin();il!=vline.end();il++) delete (*il);
	for (std::vector<TEllipse*>::iterator il=vel.begin();il!=vel.end();il++) delete (*il);

}

void HoughCartesian::draw(DCHistogramHandler* h,std::vector< std::pair<double,double> > *maxval)
{
if (CanvasHough==NULL)
	{
		CanvasHough=new TCanvas("CanvasHough","hough",800,900);
		CanvasHough->Modified();
		CanvasHough->Draw();
		CanvasHough->Divide(1,2);
		TVirtualPad* pd=CanvasHough->cd(2);
		pd->Divide(2,1);
	}
	CanvasHough->cd();
	std::stringstream htname;
	htname<<"HoughCartesian"<<theNbinTheta_<<"_"<<theNbinR_;
	TH2F* hhtx = (TH2F*) h->GetTH2(htname.str());
	TH2F* hx = (TH2F*) h->GetTH2("LocalInverse");
	TH2F* hox = (TH2F*) h->GetTH2("LocalImage");

	if (hhtx==NULL)
	{
		hhtx =(TH2F*)h->BookTH2(htname.str(),theNbinTheta_,theThetaMin_,theThetaMax_,theNbinR_,theRMin_,theRMax_);
	}
	else
		hhtx->Reset();

	if (hx==NULL)
	{
		hx =(TH2F*)h->BookTH2("LocalInverse",100,0.,0.05,100,0.,0.05);
		hox =(TH2F*)h->BookTH2("LocalImage",400,-50.,150.,200,10.,110);
		hx->SetMarkerColor(4);
	hx->SetMarkerSize(.2);
	hx->SetMarkerStyle(25);
	hox->SetMarkerColor(4);
	hox->SetMarkerSize(.2);
	hox->SetMarkerStyle(25);
	}
	else
	{
		hx->Reset();
		hox->Reset();
	}
for (int ip=0;ip<theX_.size();ip++)
		{
			double xp=theX_[ip]/(theX_[ip]*theX_[ip]+theY_[ip]*theY_[ip]);
			double yp=theY_[ip]/(theX_[ip]*theX_[ip]+theY_[ip]*theY_[ip]);
			hx->Fill(xp,yp);
			hox->Fill(theX_[ip],theY_[ip]);
		
		}
		
	
	for (uint32_t i=0;i<theNbinTheta_;i++)
	for (uint32_t j=0;j<theNbinR_;j++)
	{
		hhtx->SetBinContent(i+1,j+1,theHoughImage_[i][j]*1.);

		
	}
	
	CanvasHough->cd(1);
	hhtx->Draw("COLZ");
	//hw->Draw();
	TVirtualPad* pd=CanvasHough->cd(2);
	pd->cd(1);
	

	hx->Draw("p");
	std::vector<TLine*> vline;vline.clear();
	if (maxval!=NULL)
	{
		for (std::vector < std::pair<double,double> >::iterator ihb=maxval->begin();ihb<maxval->end();ihb++)
		{
		double theta=(*ihb).first;
		double r=(*ihb).second;
		
		double a=-1./tan(theta);
		double b=r/sin(theta);
		
		a=theta;
		b=r;
		//
		  double R=TMath::Sqrt((a*a+1.)/b/b)/2.;
		double xi=-a/2./b;
		double yi=1./2./b;
		printf("%f %f %f %f Rayon R= %f \n => %f GeV/c  Phi0 %f \n",r,theta,a,b,R,0.3*3.8*R/100.,atan(a));

		TLine* l =new TLine(0.,b,0.05,a*0.05+b);
			l->SetLineColor(2);
			l->Draw("SAME");
				vline.push_back(l);
			if (vline.size()>5) break;
		}
	}
	pd->cd(2);
	hox->Draw("p");
	std::vector<TEllipse*> vel;vel.clear();
	if (maxval!=NULL)
	{
		
		for (std::vector < std::pair<double,double> >::iterator ihb=maxval->begin();ihb<maxval->end();ihb++)
		{
		double theta=(*ihb).first;
		double r=(*ihb).second;
		
		double a=-1./tan(theta);
		double b=r/sin(theta);
		
		a=theta;
		b=r;
		//
		  double R=TMath::Sqrt((a*a+1.)/b/b)/2.;
		//
		//double R=1./2./TMath::Abs(r);
		double xi=-a/2./b;
		double yi=1./2./b;
		//printf(" %f %f %f %f Rayon R= %f => %f GeV/c\n",r,theta,a,b,R,0.3*3.8*R/100.);
		TEllipse* l =new TEllipse(xi,yi,R,R);
			l->SetLineColor(3);
			l->SetFillStyle(0);
			l->Draw("SAME");
				vel.push_back(l);
			if (vel.size()>5) break;
		}
	}
	CanvasHough->Modified();
	CanvasHough->Draw();
	//CanvasHough->WaitPrimitive();

	CanvasHough->Update();

	char c;c=getchar();putchar(c); if (c=='.') exit(0);
	for (std::vector<TLine*>::iterator il=vline.begin();il!=vline.end();il++) delete (*il);
	for (std::vector<TEllipse*>::iterator il=vel.begin();il!=vel.end();il++) delete (*il);

}
void HoughCartesian::PrintConvert(double theta,double r)
{
	double a=-1./tan(theta);
	double b=r/sin(theta);
		
		
		//
		double R=1./2./TMath::Abs(r);
		double xi=-a/2./b;
		double yi=1./2./b;
		printf(" From r=%f theta=%f a=%f b=%f  R= %f  => Pt=%f GeV/c  Phi0=%f \n",r,theta,a,b,R,0.3*3.8*R/100.,atan(a));

}


