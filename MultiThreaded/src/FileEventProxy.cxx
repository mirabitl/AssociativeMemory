#include "FileEventProxy.h"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <sys/dir.h>  
#include <sys/param.h>  
#include <stdio.h>  
#include <stdlib.h>  
#include <unistd.h>  
#include <string.h>

#include <fcntl.h>
#include <sstream>
extern  int alphasort(); //Inbuilt sorting function  
#define FALSE 0  
#define TRUE !FALSE      
int file_select_1(const struct direct *entry)  
{  
  if ((strcmp(entry->d_name, ".") == 0) || (strcmp(entry->d_name, "..") == 0))  
    return (FALSE);  
  else  
    return (TRUE);  
}

FileEventProxy::FileEventProxy(std::string dir) : theDirectory_(dir) 
{
  this->CreateDirectories();
}


void FileEventProxy::CreateDirectories()
{
  std::stringstream ss;
  ss<<"mkdir -p "<<theDirectory_;
  system(ss.str().c_str());
  ss<<"/closed";
  system(ss.str().c_str());

}

void FileEventProxy::List(std::vector<std::string> &names)
{
  names.clear();
  int count,i;  
  struct direct **files;    
  std::stringstream ss;
  ss<<theDirectory_<<"/closed";
  count = scandir(ss.str().c_str(), &files, file_select_1, alphasort);  
	
  /* If no files found, make a non-selectable menu item */  
  for (i=1; i<count+1; ++i)
    {
      std::string s(files[i-1]->d_name);
      names.push_back(s);
    }
}
void FileEventProxy::Write(std::string name,char* cbuf,uint32_t size_buf)
{
  std::stringstream sn;
  sn<<theDirectory_<<"/"<<name;
  std::stringstream snc;
  snc<<theDirectory_<<"/closed/"<<name;
  

  int fd= ::open(sn.str().c_str(),O_CREAT| O_RDWR | O_NONBLOCK,S_IRWXU);
  if (fd<0)
    {
      printf("%s No way to store to file %s:",__PRETTY_FUNCTION__,sn.str().c_str());
      //std::cout<<" No way to store to file"<<std::endl;
      return;
    }
   int ier=write(fd,cbuf,size_buf);
   if (ier!=size_buf) 
     {
       printf("%s No way to write to file %s: error %d",__PRETTY_FUNCTION__,sn.str().c_str(),ier);

       return;
     }
   ::close(fd);
  
   fd= ::open(snc.str().c_str(),O_CREAT| O_RDWR | O_NONBLOCK,S_IRWXU);
   if (fd<0)
    {
      printf("%s No way to store to file %s:",__PRETTY_FUNCTION__,snc.str().c_str());
      //std::cout<<" No way to store to file"<<std::endl;
      return;
    }
   //std::cout<<st.str().c_str()<<" "<<fd<<std::endl;
   //write(fd,b,1);
   ::close(fd);
}
void FileEventProxy::Read(std::string name,char* cbuf,uint32_t& size_buf)
{
   std::stringstream sn;
   sn<<theDirectory_<<"/"<<name;
   int fd=::open(sn.str().c_str(),O_RDONLY);
   if (fd<0) 
     {
       printf("%s No way to open to file %s: ",__PRETTY_FUNCTION__,sn.str().c_str());
       return;
     }
   size_buf=::read(fd,cbuf,0x100000);
      //printf("%d bytes read %x %d \n",size_buf,cbuf[0],cbuf[1]);
   ::close(fd);
}
void FileEventProxy::Erase(std::string name)
{
  std::stringstream sn;
  sn<<theDirectory_<<"/"<<name;
  std::stringstream snc;
  snc<<theDirectory_<<"/closed/"<<name;
  ::unlink(snc.str().c_str());
  ::unlink(sn.str().c_str());
}
