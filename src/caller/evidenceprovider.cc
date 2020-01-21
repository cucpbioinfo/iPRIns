#include "evidenceprovider.h"
#include <dirent.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

EvidenceProvider::EvidenceProvider(FileManager *fp)
{
    filepath = fp;
    calculateSizeEvidenceFromFile();
    loadEvidenceFromFile();
}

int EvidenceProvider::getEvidenceSize()
{
    // return sizeEvidenceFile;
    return evidencelist.size();
}

void EvidenceProvider::calculateSizeEvidenceFromFile()
{
    sizeEvidenceFile= 0;
    std::string line;
    std::ifstream myfile(filepath->getAllEvidencePath());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            sizeEvidenceFile++;
        }
        myfile.close();
    }
    else
        std::cout << "Unable to open file";
}

void EvidenceProvider::loadEvidenceFromFile()
{
    std::string line;
    std::ifstream myfile(filepath->getAllEvidencePath());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            Evidence e;
            e.setEvidenceByString(line);
            evidencelist.push(e);
        }
        myfile.close();
    }
    else
        std::cout << "Unable to open file";
}

Evidence EvidenceProvider::getEvidence()
{
    nextEvidence = evidencelist.front();
    evidencelist.pop();
    return nextEvidence;
}

bool EvidenceProvider::isEmpty()
{
    if (getEvidenceSize() <= 0)
    {
        return true;
    }

    return false;
}
