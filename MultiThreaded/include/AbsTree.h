#ifndef ABSTREE_H
#define ABSTREE_H
#include <string>
#include <vector>
#include <iostream>
#include "DCUtils.h"
class AbsTreeNode
{
 public:
//  friend std::ostream& operator<< (std::ostream& o, const AbsTreeNode& fred);
  AbsTreeNode(std::string name){nodeName_=name;nodeChild_.clear();parentNode_=NULL;}
  inline void addChild(AbsTreeNode *a){nodeChild_.push_back(a); a->setParentNode(this);}
  inline std::string getName()  {return nodeName_;}
  inline std::vector<AbsTreeNode*>&  getChilds(){return nodeChild_;}
  inline void setParentNode(AbsTreeNode* t){parentNode_=t;}
  inline AbsTreeNode* getParentNode(){ return parentNode_;}
  inline std::string getFullName()
  {
    
    std::string fullName=nodeName_;
    AbsTreeNode* a =this->getParentNode();
    while (a!=NULL)
      {
	fullName = a->getName()+"/"+fullName;
	a = a->getParentNode();
      }
    fullName="/"+fullName;
    return fullName;
  }
  inline void addFullName(std::string n)
    {
      //std::cout<<n<<" to be in "<<std::endl;
      std::vector<std::string> vl;
      DCUtils::String2Vector(n,vl,"/");
      AbsTreeNode* t=this;
      //printf("vl size %d \n",vl.size());
      for (unsigned int i=1;i<vl.size();i++)
	{
	  //Loop on child name
	  AbsTreeNode* son=NULL;
	   for (std::vector<AbsTreeNode*>::iterator it=t->getChilds().begin();it!=t->getChilds().end();it++)
	    {
	      //std::cout<<" Child found "<<(*it)->getName()<<std::endl;
	      if ((*it)->getName().compare(vl[i])==0) son=(*it);
	    }
	   //std::cout<<vl[i]<<" "<<(int) son<<std::endl;
	   if (son==0) 
	     {AbsTreeNode* nat=new AbsTreeNode(vl[i]);t->addChild(nat);
	       //std::cout<<"Adding "<<vl[i]<<" to "<<t->getFullName()<<std::endl;
	       t=nat;
	     }
	   else
	     t=son;
	     
	}
    }
  void toHTML(std::ostream&o)
    {
     
      o<<this->getName()<<std::endl;
      if (this->getChilds().size()>0)
	{
	  o<<"<ul>"<<std::endl;
	  for (std::vector<AbsTreeNode*>::iterator it=this->getChilds().begin();it!=this->getChilds().end();it++)
	    {

	 
	      o<<"\t  <li> "; (*it)->toHTML(o);o<<std::endl;
	    

	    
	    }
	  o<<"</ul>"<<std::endl;
	}

    }

  AbsTreeNode* findFullName(std::string name)
    {
      
      if (this->getChilds().size()>0)
	{
	  for (std::vector<AbsTreeNode*>::iterator it=this->getChilds().begin();it!=this->getChilds().end();it++)
	    {
	      
	      
	      AbsTreeNode* t= (*it)->findFullName(name);
	      if (t!=0) return t;
	    }
	}

      if (this->getFullName().compare(name)==0) 
	return this;
      else
	return 0;

    }
    



 private:
  std::string nodeName_;
  AbsTreeNode* parentNode_;
  std::vector<AbsTreeNode*> nodeChild_;
};
#ifdef AQUEJONI
  std::ostream& operator<< (std::ostream& o,  AbsTreeNode& fred)
    {
      
      AbsTreeNode* p=fred.getParentNode();
      while( p!=NULL){ o<<"     ";p=p->getParentNode();}


      if (fred.getChilds().size()>0)
	{
	  o<<"+"<<fred.getName()<<std::endl;
	 

	  for (std::vector<AbsTreeNode*>::iterator it=fred.getChilds().begin();it!=fred.getChilds().end();it++)
	    {

	 
	      o<<(**it);
	    
	    }

	
	}
      else
	  o<<"-"<<fred.getName()<<std::endl;

      return o;
    }

#endif
#endif
