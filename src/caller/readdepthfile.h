#ifndef READDEPTHFILE_H
#define READDEPTHFILE_H
#include <string>
#include "filemanager.h"
#include "readdepthhelper.h"

class ReadDepthFile
{
private:
  std::string filepath;
  FileManager *filemanager;
  // int32_t avgReadDepthFocus=0;
  std::map<int32_t, ReadDepthHelper::ReadDepthVector> mapReadDepthLineSegment;

public:
  std::string getFilePath();
    ReadDepthFile();
  void setFileManager(FileManager *filemanager);
  void loadDataToCache(std::string filepath);
  ReadDepthHelper::ReadDepthVector getReadDepthAt(int32_t number);
  std::vector<std::string> split(const std::string &s, char delimiter);
  ReadDepthHelper::ReadDepthVector findReadDepthWithFile(int32_t number,std::string filepath,int32_t scope);
  void addToMapReadDepthLineSegment(std::vector<std::string> *lineBuffer);
};

#endif