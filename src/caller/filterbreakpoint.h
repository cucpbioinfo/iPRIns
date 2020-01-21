#ifndef FilterBreakpoint_H
#define FilterBreakpoint_H
#include "evidence.h"
#include "filemanager.h"
#include "readdepthstat.h"

class FilterBreakpoint
{
private:
    FileManager *filemanager;
    int32_t vcfIdNumber = 0;
    int32_t roundConfig = 250;
    ReadDepthStat readDepthStat;
    SampleStat *samplestat;
    int minimumdivide = 4;

public:
    FilterBreakpoint();

    void execute();
    std::vector<std::string> getPathVCFFiles();
    void setFileManager(FileManager *filemanager);
    void setSampleStat(SampleStat *samplestat);
    std::vector<Evidence> getEvidenceByFilepath(std::string filepaht);
    std::vector<Evidence> getResultWithOutOverlapped(std::vector<Evidence> *master, std::vector<Evidence> *slave);
    std::vector<Evidence> getResultRemoveOverlapped(std::vector<Evidence> *master, std::vector<Evidence> *slave);
    std::vector<Evidence> getRefineResultInsertion(std::vector<Evidence> *master);

    void writeFile(std::vector<Evidence> *master);
    int32_t roundNumber(int32_t number,int32_t round);
    int32_t nextNumber(int32_t number,int32_t round);
    int32_t previousNumber(int32_t number,int32_t round);
    int getDivider(int value,int top,int down,int minimum);

    bool checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped);
};

#endif