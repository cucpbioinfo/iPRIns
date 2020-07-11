#include "readdepthanalysis.h"
#include <fstream>
#include <iostream>

ReadDepthAnalysis::ReadDepthAnalysis(FileManager *filemanager)
{
    ReadDepthAnalysis::filemanager = filemanager;
    // readDepthStat = readDepthStat();
    readDepthStat.setFilePath(filemanager);
    readDepthStat.execute();
}

void ReadDepthAnalysis::setSampleStat(SampleStat *samplestat)
{
}

int32_t ReadDepthAnalysis::getRound(int32_t x, int32_t max)
{
    return (x / max) * max;
}

std::vector<std::int32_t> ReadDepthAnalysis::getVectorRange(int32_t pos, int32_t end)
{
    std::vector<std::int32_t> temp;
    for (int32_t i = pos; i <= end; i += configRound)
    {
        temp.push_back(getRound(i, configRound));
    }
    return temp;
}

void ReadDepthAnalysis::setFocusReadDepth(int32_t pos, int32_t end, std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth)
{
    ReadDepthHelper::ReadDepthVector temprdvector;

    std::vector<std::int32_t> listrange = getVectorRange(pos, end);

    // std::vector<ReadDepthHelper::ReadDepthVector> vFocus;
    for (auto n : listrange)
    {
        ReadDepthHelper::ReadDepthVector datamodel;
        auto data = mapReadDepthLineSegment[n];
        datamodel.pos = data.pos;
        datamodel.depth = data.depth;

        focusReadDepth->push_back(datamodel);
    }
}

int ReadDepthAnalysis::getAvgReadDepth()
{
    return avgReadDepth;
}


bool ReadDepthAnalysis::filterInsertion(Evidence e)
{

    // std::cout << "getReadDepthAverageFocusArea " << getReadDepthAverageFocusArea(&startFocusReadDepth) << " " << (readDepthStat.getReadDepthByChr(e.getChr()) * 2) << std::endl;

    // if (getReadDepthAverageFocusArea(&startFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
    // {
    //     return false;
    // }

    if ((e.getMark() == "MATEUNMAPPED"))
    {
        if (e.getFrequency() <= 4)
        {
            return false;
        }

        if (e.getMaxMapQ() < 40)
        {
            return false;
        }
    }
    else if ((e.getMark() == "SINS"))
    {
    }
    else if ((e.getMark() == "SR"))
    {
    }
    else
    {
//        if (e.getFrequency() <= 1)
//        {
//            return false;
//        }

//        if (e.getMaxMapQ() < 10)
//        {
//            return false;
//        }
    }

//    if (sumStartSCL <= 3 && sumStartSCF <= 3)
//    {
//        return false;
//    }

    return true;
}


bool ReadDepthAnalysis::analyzeByEvidence(Evidence e)
{
    if (e.getMark() == "SR")
    {
        return true;
    }

    if (e.getMark() == "SDEL")
    {
        return true;
    }

    if (cachechr != e.getChr())
    {
        std::cout << e.convertToVcfString() << std::endl;

        loadDataToCache(filemanager->getReadDepthPath() + "/" + e.getChr() + ".txt");
        std::cout << "✓ : " << e.getChr() << std::endl;
        cachechr = e.getChr();
    }

    startFocusReadDepth.clear();
    endFocusReadDepth.clear();


    if (e.getVariantType() == "INS")
    {
        setFocusReadDepth(e.getPos() + e.getCiPosLeft() - configRound, e.getPos() + e.getCiPosRight() + configRound, &startFocusReadDepth);
        collectNewData();

        if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 400)
        {
            return false;
        }

        return filterInsertion(e);
    }

    return false;
}

void ReadDepthAnalysis::collectNewData()
{

    sumStartSCF = 0;
    sumStartSCL = 0;

    sumEndSCF = 0;
    sumEndSCL = 0;


    for (auto n : startFocusReadDepth)
    {

        sumStartSCF += n.SCF;
        sumStartSCL += n.SCL;

    }

    for (auto n : endFocusReadDepth)
    {

        sumEndSCF += n.SCF;
        sumEndSCL += n.SCL;

    }
}

int ReadDepthAnalysis::getSCFFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth)
{
    int count = 0;
    for (auto n : *focusReadDepth)
    {
        count += n.SCF;
    }

    return count;
}

int ReadDepthAnalysis::getSCLFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth)
{
    int count = 0;
    for (auto n : *focusReadDepth)
    {
        count += n.SCL;
    }

    return count;
}

int ReadDepthAnalysis::getReadDepthAverageFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth)
{
    int sumRD = 0;
    int count = 0;
    for (auto n : *focusReadDepth)
    {
        sumRD += n.depth;
        count++;
    }
    int avgReadDepthFocus = 0;
    if (count == 0)
    {
        avgReadDepthFocus = 0;
    }
    else
    {
        avgReadDepthFocus = int(sumRD / count);
    }
    return avgReadDepthFocus;
}

void ReadDepthAnalysis::loadDataToCache(std::string filepath)
{
    mapReadDepthLineSegment.clear();

    std::string line;
    int sumRD = 0;
    int count = 0;
    std::ifstream myfile(filepath);
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {

            ReadDepthHelper::ReadDepthVector temp;
            std::vector<std::string> token = split(line, '\t');
            temp.pos = std::stol(token.at(0), nullptr, 0);
            temp.depth = std::stoi(token.at(1));

            temp.SCF = std::stoi(token.at(2));
            temp.SCL = std::stoi(token.at(3));

            sumRD += temp.depth;
            count++;
            mapReadDepthLineSegment[temp.pos] = temp;

            // cacheReadDepthFile.push_back(temp);
        }
        myfile.close();
    }
    else
        std::cout << "Unable to open file";

    if (sumRD == 0)
    {
        avgReadDepthFocus = 0;
    }
    else
    {
        avgReadDepthFocus = int(sumRD / count);
    }

}

std::vector<std::string> ReadDepthAnalysis::split(const std::string &s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

int ReadDepthAnalysis::getDivider(int value, int top, int down, int minimum)
{
    auto returnvalue = (int)(float(value) * (float(top) / float(down)));

    if (returnvalue > minimum)
    {
        return returnvalue;
    }

    return minimum;
}