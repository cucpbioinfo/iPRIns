#include "fasta.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

Fasta::Fasta()
{
}

void Fasta::setFilePath(std::string t_filepath)
{
    filepath = t_filepath;
}

std::string Fasta::getFilePath()
{
    return filepath;
}

std::string Fasta::getSeqByOffset(int64_t start, int64_t end)
{
    std::ifstream file(filepath);
    file.seekg(start);

    int64_t sumOffset = start;
    int64_t currentPosition = file.tellg();

    std::string output;
    std::string line;
    while (std::getline(file, line))
    {
        sumOffset += line.size() + 1;

        if (sumOffset == end)
        {
            output += line.c_str();
            break;
        }
        else if (sumOffset > end)
        {
            int64_t diff = sumOffset - end;
            output += line.substr(0, line.size() - diff + 1);
            break;
        }

        output += line.c_str();
    }

    return output;
}