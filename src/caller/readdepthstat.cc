#include "readdepthstat.h"

ReadDepthStat::ReadDepthStat(FileManager *filepath)
{
    ReadDepthStat::filepath = filepath;
    
}

ReadDepthStat::ReadDepthStat() {

}

void ReadDepthStat::setFilePath(FileManager *filepath)
{
    ReadDepthStat::filepath = filepath;
}

void ReadDepthStat::execute()
{

    std::string line;
    std::ifstream myfile(filepath->getReadDepthStatPath());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            std::vector<std::string> token = split(line, '=');
            rdmap[token.at(0)] = std::stoi(token.at(1));
        }
    }
}


std::vector<std::string> ReadDepthStat::split(const std::string &s, char delimiter)
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

int ReadDepthStat::getReadDepthByChr(std::string chr) {
    return rdmap[chr];
}
