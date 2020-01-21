#ifndef READDEPTHSTAT_H
#define READDEPTHSTAT_H

#include <string>
#include <vector>
#include "readdepthhelper.h"
#include "filemanager.h"
#include <fstream>

class ReadDepthStat {
private:
    FileManager *filepath;
    std::map<std::string, int> rdmap;
    
public:
    ReadDepthStat(FileManager *filepath);
    ReadDepthStat();
    void execute();
    void setFilePath(FileManager *filepath);
    std::vector<std::string> split(const std::string &s, char delimiter);
    int getReadDepthByChr(std::string chr);
};

#endif