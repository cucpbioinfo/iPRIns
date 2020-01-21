#ifndef READDEPTHHELPER_H
#define READDEPTHHELPER_H
#include <cstdint>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "filemanager.h"
#include "evidencefinder.h"

class ReadDepthHelper
{
private:
    FileManager *filepath;

    int avgReadDepth = 0;

public:
    int getAvgReadDepth();

    void setAvgReadDepth(int avgReadDepth);
    struct ReadDepthVector
    {
        int32_t pos = 0;
        int depth = 0;

        int SCF = 0;
        int SCL = 0;


        bool operator<(const ReadDepthVector &rhs) const
        {
            return (pos < rhs.pos);
        }
    };

private:
    std::string outputpath;
    std::string target_chromosome;
    int avgReaddepth = 0;
    int32_t range = 0;

    struct VariantRangeRD
    {
        int32_t pos = 0;
        int32_t end = 0;

        bool operator<(const ReadDepthVector &rhs) const
        {
            return (pos < rhs.pos);
        }
    };

    std::vector<ReadDepthVector> vecRDLine;

    std::vector<VariantRangeRD> LowRD;

    std::vector<std::string> split(const std::string &s, char delimiter);

protected:
public:
    ReadDepthHelper(FileManager *filepath);
    ReadDepthHelper();
    void setReadDepthMap(std::map<int32_t, EvidenceFinder::ReadDepthDetail> rdmap);
    void writeReadDepthLineFile(std::string path);
    void findEvidence();
    void setRange(int32_t range);
    void calculateAvgReaddepth();
    void setTargetChromosome(std::string chr);
    void setOutputPath(std::string path);
    void loadReadDepthFile(std::string path);
    std::vector<std::string> getVariantVcfFormat();
    void writeVcf(std::string path);
    bool isRangeDisorderByMorethanRD(int32_t pos, int32_t end, int readdepth);
    void writeReadDepthStat(std::string path);
};

#endif