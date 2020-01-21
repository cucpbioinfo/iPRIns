#ifndef SPLITREAD_H
#define SPLITREAD_H
#include <string>
#include "readparser.h"
#include "samplestat.h"
#include "refiningsv.h"
#include "evidence.h"
#include "filemanager.h"

class SplitRead
{
private:
  ReadParser *readparser;
  SampleStat *samplestat;
  std::string chrname;
  std::vector<ReadParser::SATag> satag;
  FileManager *filepath;


  std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> mapSmallINS;
  std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> mapINS;


  int vcfIdNumber = 0;

public:
  SplitRead(std::string chrname,ReadParser *readparser, SampleStat *samplestate,FileManager *filepath);
  void findInsertionInRead();
  void updateRead();
  void printResult();
  void removeDuplicateResult(std::vector<Evidence> *vec);
  bool checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped);
  int writeFile(Evidence vr);

  void printSmallInsertion();
  void findInsertion();

  void printInsertion();
  std::vector<Evidence> convertMapToEvidenceList(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *mapSV, std::string svtype,std::string mark);
  void mergeEvidence(std::vector<Evidence> *elist);
  void setAllCIEvidence(std::vector<Evidence> *elist,int32_t rangePos);
  void filterEvidenceList(std::vector<Evidence> *elist);
  void filterLengthMinEvidenceList(std::vector<Evidence> *elist,int32_t min);
  void filterLengthMaxEvidenceList(std::vector<Evidence> *elist,int32_t min);
    void filterMapQLowerThan(uint8_t mapq,std::vector<Evidence> *elist);
     void filterFrequencyLowerThan(int number,std::vector<Evidence> *elist);
    int getDivider(int value, int top, int down, int minimum);
    std::vector<Evidence> convertMapToEvidenceListTRA(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *mapSV, std::string svtype, std::string mark);

  
};

#endif