#ifndef EVIDENCEPROVIDER_H
#define EVIDENCEPROVIDER_H
#include "evidence.h"
#include "filemanager.h"
#include <queue>

class EvidenceProvider
{

  private:
    FileManager *filepath;

    Evidence nextEvidence;
    void loadEvidenceFromFile();
    std::queue<Evidence> evidencelist;
     void calculateSizeEvidenceFromFile();

     int sizeEvidenceFile=0;
  public:
    EvidenceProvider(FileManager *filepath);
    int getEvidenceSize();
   
    Evidence getEvidence();
    bool isEmpty();
};

#endif