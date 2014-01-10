#ifndef _FileEventProxy_H_
#define _FileEventProxy_H_
#include <string>
#include <vector>
#include <stdint.h>

class FileEventProxy
{
 public:
  FileEventProxy(std::string directory="/dev/shm/Events");
  void CreateDirectories();
  void List(std::vector<std::string> &names,std::string pattern="*");
  void Write(std::string name,char* buf,uint32_t size);
  void Read(std::string name,char* buf,uint32_t& size);
  void Erase(std::string name);
 protected:
  std::string theDirectory_;

};
#endif
