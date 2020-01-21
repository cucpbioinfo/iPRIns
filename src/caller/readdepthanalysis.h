#ifndef READDEPTHANALYSIS_H
#define READDEPTHANALYSIS_H

#include <string>
#include <vector>
#include "filemanager.h"
#include "readdepthhelper.h"
#include <sstream>
#include "evidence.h"
#include "samplestat.h"
#include <map>
#include "readdepthstat.h"

class ReadDepthAnalysis {
private:
    FileManager *filemanager;
    std::string cachechr;
    int avgReadDepthFocus=0;
    int avgReadDepth = 0;
    std::map<int32_t, ReadDepthHelper::ReadDepthVector> mapReadDepthLineSegment;
    // std::vector<ReadDepthHelper::ReadDepthVector> cacheReadDepthFile;
    std::vector<ReadDepthHelper::ReadDepthVector> startFocusReadDepth;
    std::vector<ReadDepthHelper::ReadDepthVector> endFocusReadDepth;
    std::vector<std::string> split(const std::string &s, char delimiter);
    SampleStat *samplestat;
    uint32_t configRound = 250;

    ReadDepthStat readDepthStat;
    int minimumdivide = 2;
    int getDivider(int value,int top,int down,int minimum);

private:

    int sumStartSCF = 0;
    int sumStartSCL = 0;


    int sumEndSCF = 0;
    int sumEndSCL = 0;



public:
    void collectNewData();
    ReadDepthAnalysis(FileManager *filemanager);
    int32_t getRound(int32_t x, int32_t max);
    std::vector<std::int32_t> getVectorRange(int32_t pos, int32_t end);
    bool analyzeByEvidence(Evidence e);
    void setMedianReadDepth(int readdepth);
    int getReadDepthAverageFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth);
    void setSampleStat(SampleStat *samplestat);
    void loadDataToCache(std::string filepath);
    void loadEvidenceFromFile();
    int getAvgReadDepth();
    void setFocusReadDepth(int32_t pos, int32_t end,std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth);
    int getNumberReadDepthVector(std::vector<ReadDepthHelper::ReadDepthVector> focus,int positionnumber);
    bool filterInsertion(Evidence e);
    void loadavgReadDepthFocusStat();


    int getSCFFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth);
    int getSCLFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth);
};

#endif