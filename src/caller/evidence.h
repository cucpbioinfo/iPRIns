#ifndef EVIDENCE_H
#define EVIDENCE_H

#include <string>
#include <stdint.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

class Evidence
{
private:
  std::string variantType;
  std::string chr;
  std::string mchr;
  int32_t posDiscordantRead = 0;
  int32_t endDiscordantRead = 0;
  int32_t lastposDiscordantRead = 0;
  int32_t lastendDiscordantRead = 0;

  int32_t pos = 0;
  int32_t end = 0;
  int32_t lastpos = 0;
  int32_t lastend = 0;

  int32_t ciPosLeft = 0;
  int32_t ciPosRight = 0;

  bool foundEvidenceAtStart = false;
  bool foundEvidenceAtEnd = false;

  std::vector<std::string> multipleEndChromosome;

  std::string mark;
  std::vector<uint8_t> rpmapqlist;

public:
  void setRPMapQ(std::vector<uint8_t> rpmapq);
  std::vector<uint8_t> *getRPMapQ();

  // void setMark(std::string mark);
  // std::string getMark();
  int LNGMATCH = 0;

  void setMark(std::string mark);
  std::string getMark();
  std::vector<std::string> &getMultipleEndChromosome();

  void setMultipleEndChromosome(const std::vector<std::string> &multipleEndChromosome);
  void setMapQList(std::vector<uint8_t> mapqs);
  int getNumberOfZeroMapQ();
   int getNumberOfZeroRPMapQ();

public:
  //    bool isFoundEvidenceAtStart() const;
  //
  //    void setFoundEvidenceAtStart(bool foundEvidenceAtStart);
  //
  //    bool isFoundEvidenceAtEnd() const;
  //
  //    void setFoundEvidenceAtEnd(bool foundEvidenceAtEnd);
  std::string convertMapQlistToCommaString(std::vector<uint8_t> *mapqs);

private:
  int32_t ciEndLeft = 0;
  int32_t ciEndRight = 0;

  int frequency = 0;
  bool forwardDirection;
  bool backwardDirection;

  struct associateRead
  {
    int32_t posDiscordantRead = 0;
    int32_t mateposDiscordantRead = 0;
    // int count=0;
  };

  std::vector<associateRead> associateReadLists;

  //    int32_t calculateSDEndAssociateReadListsRead(std::vector<associateRead> *m_associateReadLists);
  //    int calculateNumberOfAbnormalSDAssociateReadListsRead(std::vector<associateRead> *m_associateReadLists,int32_t basediff);

  //    std::vector<associateRead> errorAssociateReadLists;

  //    std::string convertEndAssociateReadListsToCommaString();
  bool QuailtyPass = false;
  std::string ID;
  std::string ref;
  std::string alt;
  std::string info;
  uint8_t qual = 0;
  std::string filter;
  std::vector<uint8_t> mapqlist;
  std::string comment;

public:
  int NumberforwardDirection;
  int NumberbackwardDirection;

  void setID(std::string id);
  int32_t getCiPosLeft();

  void setCiPosLeft(int32_t ciPosLeft);

  int32_t getCiPosRight();

  void setCiPosRight(int32_t ciPosRight);

  int32_t getCiEndLeft();

  void setCiEndLeft(int32_t ciEndLeft);

  int32_t getCiEndRight();

  void setCiEndRight(int32_t ciEndRight);
  void setQuailtyPass(bool QuailtyPass);
  bool isQuailtyPass() const;
  std::string getResultVcfFormatString();
  std::string getID();
  std::string getRef();
  std::string getAlt();
  uint8_t getQual();

  std::string getFilter();

  std::string getInfoString();
  // std::string getComment();
  // void setComment(std::string comment);

public:
  Evidence();

  ~Evidence();

  int getSizeErrorAssociateReadLists();

  std::string getVariantType();

  void setVariantType(std::string varianttype);

  void setPosDiscordantRead(int32_t pos);

  void setEndDiscordantRead(int32_t end);

  void setLastPosDiscordantRead(int32_t pos);

  void setLastEndDiscordantRead(int32_t end);

  void setFrequency(int count);

  void incrementFrequency();

  void setChr(std::string chr);

  void setEndChr(std::string matechr);

  void setForwardDirection(bool isforward);

  void addMapQ(uint8_t mapq);

  void addAssociateRead(int32_t pos, int32_t matepos);
  //    void addErrorAssociateRead(int32_t pos, int32_t matepos);

  std::vector<uint8_t> *getMapQVector();

  int32_t getPosDiscordantRead();

  int32_t getEndDiscordantRead();

  int getFrequency();

  bool getForwardDirection();

  void setBackwardDirection(bool isbackward);
  bool getBackwardDirection();
  std::string getChr();

  std::string getEndChr();

  int32_t getLastPosDiscordantRead();

  int32_t getLastEndDiscordantRead();

  std::string getValuebyText(std::string txt);

  std::string getKeybyText(std::string n);

  std::string convertToVcfString();

  void setEvidenceByString(std::string line);

  bool operator<(const Evidence &otherEvidence) const
  {
    return pos < otherEvidence.pos;
  }

  int32_t getPos();

  void setPos(int32_t pos);

  int32_t getEnd();

  void setEnd(int32_t end);

  int32_t getLastPos();

  void setLastPos(int32_t lastpos);

  int32_t getLastEnd();

  void setLastEnd(int32_t lastend);

  uint8_t getMaxMapQ();
  uint8_t getMinMapQ();
  uint8_t getAvgMapQ();

  uint8_t getMaxRPMapQ();
  uint8_t getMinRPMapQ();
  uint8_t getAvgRPMapQ();
  int getNumberOfRP();

  int32_t getSvLength();
  std::string getSVType();

  bool haveSomeMapQMoreThan(uint8_t qual);
  bool haveSomeMapQLessThan(uint8_t qual);

  struct AlternativeSA
  {
    std::string chr;
    int32_t pos = 0;
  };
  std::vector<AlternativeSA> AltSA;

  void addAlterSA(std::string chr,int32_t pos);
  int getAltSASize();
  int32_t evidencefrom = 0;
  int32_t getEvidencePos();
  void setEvidenceFrom(int32_t pos);
  int countEvidencePosNear(int32_t pos,int32_t overlapped);
  int countDiffEvidencePos(int32_t overlapped);
  bool checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped);

  //    std::string getIsFoundEvidenceAtStartString();
  //    std::string getIsFoundEvidenceAtEndString();
};

#endif