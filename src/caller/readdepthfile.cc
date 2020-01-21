#include "readdepthfile.h"
#include <map>

ReadDepthFile::ReadDepthFile()
{
}


ReadDepthHelper::ReadDepthVector ReadDepthFile::getReadDepthAt(int32_t number)
{
    return mapReadDepthLineSegment[number];
}

ReadDepthHelper::ReadDepthVector ReadDepthFile::findReadDepthWithFile(int32_t number, std::string filepath, int32_t scope)
{

    std::string line;
    std::ifstream myfile(filepath);
    int count = 0;
    ReadDepthHelper::ReadDepthVector temp;
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            count += scope;
            if (count != number)
            {
                continue;
            }

            std::vector<std::string> token = split(line, '\t');
            temp.pos = std::stol(token.at(0), nullptr, 0);
            temp.depth = std::stoi(token.at(1));

            temp.SCF = std::stoi(token.at(13));
            temp.SCL = std::stoi(token.at(14));
            break;
        }
        myfile.close();
    }

    return temp;
}

std::vector<std::string> ReadDepthFile::split(const std::string &s, char delimiter)
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

std::string ReadDepthFile::getFilePath()
{
    return filepath;
}

void ReadDepthFile::loadDataToCache(std::string filepath)
{
    if (ReadDepthFile::filepath == filepath)
    {
        return;
    }
    ReadDepthFile::filepath = filepath;

    mapReadDepthLineSegment.clear();
    std::string line;
    // int32_t sumRD = 0;
    // int count = 0;

    std::cout << "loadDataToCache : " << filepath << std::endl; 
    std::ifstream myfile(filepath);
    if (myfile.is_open())
    {
        std::vector<std::string> lineBuffer;
        while (getline(myfile, line))
        {
            lineBuffer.push_back(line);

            if (line.size() > 10000)
            {
                addToMapReadDepthLineSegment(&lineBuffer);
                lineBuffer.clear();
            }

            // sumRD += temp.depth;
            // count++;
        }

        if (lineBuffer.size() != 0)
        {
            addToMapReadDepthLineSegment(&lineBuffer);
            lineBuffer.clear();
        }

        myfile.close();
    }
    else
        std::cout << "Unable to open file";

}

void ReadDepthFile::addToMapReadDepthLineSegment(std::vector<std::string> *lineBuffer)
{
    for (auto line : *lineBuffer)
    {
        ReadDepthHelper::ReadDepthVector temp;
        std::vector<std::string> token = split(line, '\t');
        temp.pos = std::stol(token.at(0), nullptr, 0);
        temp.depth = std::stoi(token.at(1));

        temp.SCF = std::stoi(token.at(13));
        temp.SCL = std::stoi(token.at(14));

        mapReadDepthLineSegment[temp.pos] = temp;
    }
}