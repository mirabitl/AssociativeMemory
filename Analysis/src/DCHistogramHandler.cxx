#include "DCHistogramHandler.h"
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <sstream>
#include "DCUtils.h"
#include "TSystem.h"
#include <iomanip>


#include "TBufferXML.h"
DCHistogramHandler* DCHistogramHandler::_me = 0 ;


DCHistogramHandler* DCHistogramHandler::instance() {
  
  if( _me == 0 ) 
    _me = new DCHistogramHandler ;
  
  return _me ;
}  


DCHistogramHandler::DCHistogramHandler()
{
  mapH1.clear();
  mapH2.clear();
  top = new AbsTreeNode("Top");
  sem_init(&theMutex_,0,1);
 }

void DCHistogramHandler::Lock()
{
  sem_wait(&theMutex_);
}
void DCHistogramHandler::UnLock()
{
  sem_post(&theMutex_);
}

TH1* DCHistogramHandler::BookTH1(std::string name,int nbinx,double xmin,double xmax)
{
  top->addFullName(name);
  TH1F* h =  new TH1F(name.c_str(),name.c_str(),nbinx,xmin,xmax);
  h->SetDirectory(0);
  std::pair<std::string,TH1*> pr(name,h);
  mapH1.insert(pr);
  //std::cout<<name<<" Booked"<<std::endl;
  return h;
}

TH2* DCHistogramHandler::BookTH2(std::string name,int nbinx,double xmin,double xmax,int nbiny,double ymin,double ymax)
{
  top->addFullName(name);
  TH2F* h =  new TH2F(name.c_str(),name.c_str(),nbinx,xmin,xmax,nbiny,ymin,ymax);
  h->SetDirectory(0);
  std::pair<std::string,TH2*> pr(name,h);
  mapH2.insert(pr);
  //std::cout<<name<<" Booked"<<std::endl;
  return h;
}

TH3* DCHistogramHandler::BookTH3(std::string name,int nbinx,double xmin,double xmax,int nbiny,double ymin,double ymax,int nbinz,double zmin,double zmax)
{
  top->addFullName(name);
  TH3C* h =  new TH3C(name.c_str(),name.c_str(),nbinx,xmin,xmax,nbiny,ymin,ymax,nbinz, zmin, zmax);
  h->SetDirectory(0);
  std::pair<std::string,TH3*> pr(name,h);
  mapH3.insert(pr);
  //std::cout<<name<<" Booked"<<std::endl;
  return h;
}

TProfile* DCHistogramHandler::BookProfile(std::string name,int nbinx,double xmin,double xmax,double ymin,double ymax)
{
  top->addFullName(name);
  TProfile* h;
  if (ymin<ymax) h=  new TProfile(name.c_str(),name.c_str(),nbinx,xmin,xmax,ymin,ymax,"S");
  else  h=  new TProfile(name.c_str(),name.c_str(),nbinx,xmin,xmax);
  h->SetDirectory(0);
  h->Sumw2();
  std::pair<std::string,TH1*> pr(name,h);
  //std::cout<<name<<" Booked"<<std::endl;
  mapH1.insert(pr);
  return h;
}
TProfile2D* DCHistogramHandler::BookProfile2D(std::string name,int nbinx,double xmin,double xmax,int nbiny, double ymin,double ymax, double zmin, double zmax)
{
  top->addFullName(name);
  TProfile2D* h;
  if (zmin<zmax)h=  new TProfile2D(name.c_str(),name.c_str(),nbinx,xmin,xmax,nbiny,ymin,ymax,zmin,zmax);
  else h=new TProfile2D(name.c_str(),name.c_str(),nbinx,xmin,xmax,nbiny,ymin,ymax);
  h->SetDirectory(0);
  h->Sumw2();
  std::pair<std::string,TH2*> pr(name,h);
  //std::cout<<name<<" Booked"<<std::endl;
  mapH2.insert(pr);
  return h;
}


void DCHistogramHandler::writeHistograms(std::string name)
{
  
  TFile f(name.c_str(),"RECREATE","RAW DATA");
  f.cd();
 
  for (std::map<std::string,TH1*>::iterator iter= mapH1.begin();iter!=mapH1.end();iter++)
    {
      f.cd();
     //std::cout<<" writing "<<iter->first<<std::endl;
     std::vector<std::string> vdir;
     DCUtils::String2Vector(iter->first,vdir,"/");
     std::string dirname(""),lastname("");
     if (vdir.size()>0)
       {
       for (unsigned int i=1;i<vdir.size()-1;i++)
	 {


	   dirname+="/"+vdir[i];
	   f.cd();
	   TDirectory* det=f.GetDirectory(dirname.c_str());
	   if (det==0)
	     {
	       f.cd(lastname.c_str());
	       gDirectory->cd();
	       TDirectory* detn = gDirectory->mkdir(vdir[i].c_str(),vdir[i].c_str());
	       if (detn!=0) detn->cd();
	       det=f.GetDirectory(dirname.c_str());
	     }
	   if (det!=0) det->cd();
	   lastname=dirname;
	 }
       iter->second->SetName(vdir[vdir.size()-1].c_str());       
       }


     iter->second->Write();
     delete iter->second;
    }

  for (std::map<std::string,TH2*>::iterator iter= mapH2.begin();iter!=mapH2.end();iter++)
    {

      f.cd();
     //std::cout<<" writing "<<iter->first<<std::endl;
     std::vector<std::string> vdir;
     DCUtils::String2Vector(iter->first,vdir,"/");
     std::string dirname(""),lastname("");
     if (vdir.size()>0)
       {
       for (unsigned int i=1;i<vdir.size()-1;i++)
	 {


	   dirname+="/"+vdir[i];
	   f.cd();
	   TDirectory* det=f.GetDirectory(dirname.c_str());
	   if (det==0)
	     {
	       f.cd(lastname.c_str());
	       gDirectory->cd();
	       TDirectory* detn = gDirectory->mkdir(vdir[i].c_str(),vdir[i].c_str());
	       if (detn!=0) detn->cd();
	       det=f.GetDirectory(dirname.c_str());
	     }
	   if (det!=0) det->cd();
	   lastname=dirname;
	 }
       iter->second->SetName(vdir[vdir.size()-1].c_str());
       }    
     iter->second->Write();
     delete iter->second;
    }

  for (std::map<std::string,TH3*>::iterator iter= mapH3.begin();iter!=mapH3.end();iter++)
    {

      f.cd();
     //std::cout<<" writing "<<iter->first<<std::endl;
     std::vector<std::string> vdir;
     DCUtils::String2Vector(iter->first,vdir,"/");
     std::string dirname(""),lastname("");
     if (vdir.size()>0)
       {
       for (unsigned int i=1;i<vdir.size()-1;i++)
	 {


	   dirname+="/"+vdir[i];
	   f.cd();
	   TDirectory* det=f.GetDirectory(dirname.c_str());
	   if (det==0)
	     {
	       f.cd(lastname.c_str());
	       gDirectory->cd();
	       TDirectory* detn = gDirectory->mkdir(vdir[i].c_str(),vdir[i].c_str());
	       if (detn!=0) detn->cd();
	       det=f.GetDirectory(dirname.c_str());
	     }
	   if (det!=0) det->cd();
	   lastname=dirname;
	 }
       }







      iter->second->Write();
      delete iter->second;
    }





  mapH1.clear();
  mapH2.clear();
  mapH3.clear();
  f.Write();
  f.Close();
}
 
 
#include <sys/stat.h>
#include <sys/types.h>
#include "TBufferXML.h"
#include <iostream>
#include <fstream>

void DCHistogramHandler::writeXML(std::string  path)
{
 	std::string curpath=path;
		std::cout<<path<<std::endl;

	mode_t process_mask = umask(0);
	mkdir(curpath.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);
 	for (std::map<std::string,TH1*>::iterator iter= mapH1.begin();iter!=mapH1.end();iter++)
	{
		//std::cout<<" writing "<<iter->first<<std::endl;
		std::vector<std::string> vdir;
		DCUtils::String2Vector(iter->first,vdir,"/");
		std::string dirname(""),lastname("");
		std::string curpath=path;
		process_mask = umask(0);
		mkdir(curpath.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);	
		umask(process_mask);
		if (vdir.size()>0)
		{
			for (unsigned int i=1;i<vdir.size()-1;i++)
			{


				curpath+="/"+vdir[i];
				process_mask = umask(0);
				mkdir(curpath.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);	
				umask(process_mask);
				
			}
			lastname=curpath+"/"+vdir[vdir.size()-1];
		}
		else
			lastname=curpath+"/"+iter->first;

		//std::cout<<lastname<<std::endl;
		std::ostringstream xmlstr;
	  
	  xmlstr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
	  
	  xmlstr<<TBufferXML::ConvertToXML(iter->second);
	  ofstream myfile;
		myfile.open (lastname.c_str());
		myfile << xmlstr.str();
		myfile.close();
	  
		
		
		//iter->second->Write();
		//delete iter->second;
	}
	curpath=path;
	for (std::map<std::string,TH2*>::iterator iter= mapH2.begin();iter!=mapH2.end();iter++)
	{
		//std::cout<<" writing "<<iter->first<<std::endl;
		std::vector<std::string> vdir;
		DCUtils::String2Vector(iter->first,vdir,"/");
		std::string dirname(""),lastname("");
		std::string curpath=path;
		process_mask = umask(0);
		mkdir(curpath.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);	
		umask(process_mask);
		if (vdir.size()>0)
		{
			for (unsigned int i=1;i<vdir.size()-1;i++)
			{


				curpath+="/"+vdir[i];
				process_mask = umask(0);
				mkdir(curpath.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);	
				umask(process_mask);
				
			}
			lastname=curpath+"/"+vdir[vdir.size()-1];
		}
		else
			lastname=curpath+"/"+iter->first;

		//std::cout<<lastname<<std::endl;
		std::ostringstream xmlstr;
	  
	  xmlstr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
	  
	  xmlstr<<TBufferXML::ConvertToXML(iter->second);
	  ofstream myfile;
		myfile.open (lastname.c_str());
		myfile << xmlstr.str();
		myfile.close();
	  
		
		
		//iter->second->Write();
		//delete iter->second;
	}

 
}


void DCHistogramHandler::ListHisto(std::string reg,int type,std::vector<std::string> &found)
{
  found.clear();
  regex_t preg; int cflags=REG_EXTENDED;cflags=0;
  regcomp(&preg,reg.c_str(),cflags);

  //std::cout<<"On entre dans ListHisto "<<reg<<std::endl;
  //regmatch_t pmatch[256]; int eflags=0;size_t nmatch=0;
  if ( type ==1)
    for (std::map<std::string,TH1*>::const_iterator iter = mapH1.begin(); iter!=mapH1.end();++iter)
      {
	int ier =regexec(&preg, iter->first.c_str(),0,NULL,0);
	if (ier ==0)
	  {
	    //std::cout<<"GetTH1(\""<<iter->first<<"\")"<<std::endl;
	  found.push_back(iter->first);
	  }
      }
  if ( type ==2)
    for (std::map<std::string,TH2*>::const_iterator iter = mapH2.begin(); iter!=mapH2.end();++iter)
      {
	int ier =regexec(&preg, iter->first.c_str(),0,NULL,0);
	if (ier ==0)
	  {
	    //std::cout<<"GetTH2(\""<<iter->first<<"\")"<<std::endl;
	    found.push_back(iter->first);
	  }
      }
  if ( type ==3)
    for (std::map<std::string,TH3*>::const_iterator iter = mapH3.begin(); iter!=mapH3.end();++iter)
      {
	int ier =regexec(&preg, iter->first.c_str(),0,NULL,0);
	if (ier ==0)
	  {
	    //std::cout<<"GetTH2(\""<<iter->first<<"\")"<<std::endl;
	    found.push_back(iter->first);
	  }
      }

  //std::cout<<"On sort de ListHisto"<<std::endl;
}

TCanvas* DCHistogramHandler::DrawRegexp(std::string reg,bool same)
{

  std::vector<std::string> found;
  std::vector<TH1*> vh1;
  std::vector<TH2*> vh2;
  std::vector<TH3*> vh3;
  ListHisto(reg,1,found);
  for (unsigned int i=0;i<found.size();i++)
    vh1.push_back(GetTH1(found[i]));
  ListHisto(reg,2,found);
  for (unsigned int i=0;i<found.size();i++)
    vh2.push_back(GetTH2(found[i]));

  ListHisto(reg,3,found);
  for (unsigned int i=0;i<found.size();i++)
    vh3.push_back(GetTH3(found[i]));

  if ((vh1.size()+vh2.size()+vh3.size())==0) 
    {
      std::string serr=" No histogram found for "+reg;
      std::cout<<serr<<std::endl;
      return NULL;

    }
  //std::cout<<vh1.size()<<" TH1 Found and "<<vh2.size()<<" Th2"<<std::endl;
  std::stringstream name("Canvas");
  name<<reg;
  int ncol=TMath::Nint(TMath::Sqrt(vh1.size()+vh2.size()+1));int nrow=(vh1.size()+vh2.size())/ncol+1;
  int nhist=vh1.size()+vh2.size()+vh3.size();
  if (nhist==1)
    {
      ncol=1;
      nrow=1;
    }
  //std::cout<<"DrawRegexp :: "<<vh1.size()+vh2.size()<< " "<<ncol<<" "<<nrow<<std::endl;
  TCanvas* c1 = new TCanvas(name.str().c_str(),reg.c_str(),700,700);
  //std::cout<<"Createion du canvas faite "<<std::endl;
  if (!same)
    {
      c1->Divide(ncol,nrow);
      int ipad =1;
      gROOT->GetStyle("Default")->SetOptStat(1111);
      for (unsigned int i=0;i<vh1.size();i++)
	{
	  c1->cd(ipad);
	  if (vh1[i]== NULL) continue;
	  vh1[i]->Draw();
	  //std::cout<<"Histo plot"<<i<<std::endl;
	  ipad++;
	}
      if (vh2.size()>0)
	{
	  gROOT->GetStyle("Default")->SetOptStat(0);
	  gROOT->GetStyle("Default")->SetPalette(1,0);
	}
      for (unsigned int i=0;i<vh2.size();i++)
	{
	  c1->cd(ipad);
	  if (vh2[i]== NULL) continue;
	  vh2[i]->Draw("colz");
	  ipad++;
	}
      for (unsigned int i=0;i<vh3.size();i++)
	{
	  c1->cd(ipad);
	  if (vh3[i]== NULL) continue;
	  vh3[i]->Draw("box");
	  ipad++;
	}
    }
  else
    {
      // work only for One D histo
      bool first = true;
      gROOT->GetStyle("Default")->SetOptStat(1111);
      for (unsigned int i=0;i<vh1.size();i++)
	{
	  
	  if (vh1[i]== NULL) continue;
	  if (first)
	    {
	      vh1[i]->Draw();
	      first=false;
	    }
	  else
	    vh1[i]->Draw("SAME");
	  //std::cout<<"Histo plot"<<i<<std::endl;
	}

    }
  //std::cout<<"End of DrawRegexp"<<std::endl;
  return c1;

}
std::string DCHistogramHandler::getXMLHistoList()
{
  std::ostringstream xmlstr;
  
  xmlstr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  xmlstr << "<response success=\"true\">" << std::endl;
  for (std::map<std::string,TH1*>::const_iterator iter = mapH1.begin(); iter!=mapH1.end();++iter)
    {
      xmlstr<<"<Histo>"<<iter->first<<"</Histo>"<<std::endl;
    }
  for (std::map<std::string,TH2*>::const_iterator iter = mapH2.begin(); iter!=mapH2.end();++iter)
    {
      xmlstr<<"<Histo>"<<iter->first<<"</Histo>"<<std::endl;
    }
  for (std::map<std::string,TH3*>::const_iterator iter = mapH3.begin(); iter!=mapH3.end();++iter)
    {
      xmlstr<<"<Histo>"<<iter->first<<"</Histo>"<<std::endl;
    }
  xmlstr << "</response>" ;
  return  xmlstr.str();
  
  
}

std::string DCHistogramHandler::getXMLHisto(std::string hname)
{
  std::ostringstream xmlstr;
  //  xmlstr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  TH1* h1= this->GetTH1(hname);
  if (h1!=0)
    xmlstr<<TBufferXML::ConvertToXML(h1);
  else
    {
      TH2* h2= this->GetTH2(hname);
      if (h2!=0)
	xmlstr<<TBufferXML::ConvertToXML(h2);
      else
	{
	  TH3* h3= this->GetTH3(hname);
	  if (h3!=0)
	    xmlstr<<TBufferXML::ConvertToXML(h3);
	}
    }
 return xmlstr.str();
}
