#ifndef READPARSER_H
#define READPARSER_H
#include "samplestat.h"
#include "filemanager.h"
#include <stdint.h>
#include <iostream>
#include <string>
#include <map>
#include <vector>
/*-------------------------------------------------------
This class has the idea of having one Read pointer and then updating the pointer continuously.
this class require:
  1. bam header from htslib (bam_hdr_t)
    - this class require because use for find chromosome name
  2. pointer of read from htslib (bam1_t)

---------------------------------------------------------*/

class ReadParser
{
private:
  bam1_t *source_bamread;
  bam_hdr_t *bam_header;

  

public:
class Cigar
  {
  private:
    // int operation;
    char operation_n;
    int length = 0;

  public:
    int getLength()
    {
      return length;
    }

    char getOperatorName()
    {
      return operation_n;
    }

    void setOperatorName(char oper) {
      operation_n = oper;
    }

    void setOperatorName(int oper)
    {
      switch (oper)
      {

      case 0:
        operation_n = 'M';
        break;
      case 1:
        operation_n = 'I';
        break;
      case 2:
        operation_n = 'D';
        break;
      case 3:
        operation_n = 'N';
        break;
      case 4:
        operation_n = 'S';
        break;
      case 5:
        operation_n = 'H';
        break;
      case 6:
        operation_n = 'P';
        break;
      case 7:
        operation_n = '=';
        break;
      case 8:
        operation_n = 'X';
        break;
      default:
        std::cout << "readparser error setOperatorName() no oper" << std::endl;
        exit(0);
      }

    }
    void setLength(int length_m)
    {
      length = length_m;
    }
  };

  int getLastToStartMissMatchPosMD();
  int getStartToEndMissMatchPosMD();
  struct AlignMD
  {
    char operate;
    int size;
  };

  std::vector<ReadParser::AlignMD> getAlignMD();

protected:
public:
  ~ReadParser();
  ReadParser();
  void setBamRead(bam1_t *bam_read);
  void setBamHeader(bam_hdr_t *bam_header);
  int32_t getPos();
  int32_t getEnd();

  int32_t getMatePos();
  int32_t getMateEnd();
  int32_t getInsertSize();

  int32_t getPosOfSeq();
  int32_t getEndOfSeq();

  bool isMateUnmapped();
  bool isUnmapped();
  std::string getChromosomeNameString();
  std::string getMateChromosomeNameString();

  bool isPairOnSameChromosome();

  uint32_t getLengthCigar();

  void printCigar();

  bool hasLastCigarSoftclipped();
  bool hasFirstCigarSoftclipped();

  bool hasLastCigarHardclipped();
  bool hasFirstCigarHardclipped();

  std::string getSequence();
  int *getBaseQuality();
  uint8_t getMapQuality();
  int32_t getLengthSequence();

  std::vector<ReadParser::Cigar> getCigar();

  bool isReverse();
  bool isMateReverse();

  bool isFirstRead();
  bool isSecondRead();
  std::string getReverseComplement(std::string word);

  void replaceToReverseComplement(std::string *nucs);
  std::string getEdgeSeqFromStartSeq(uint32_t number);
  std::string getEdgeSeqFromEndSeq(uint32_t number);

  bool isNotPassingFilters();
  bool isPCR();
  bool isSupplementaryAlignment();

  int getSecondLastToStartMissMatchPosMD();
  int getSecondStartToEndMissMatchPosMD();

  std::string getSoftClippedSequenceStart();
  std::string getSoftClippedSequenceEnd();

  struct SATag
  {
    std::string chrname;
    int32_t pos=0;
    std::string strand;
    std::vector<ReadParser::Cigar> cigar;
    uint8_t mapQ=0;
    int NM=0;
  };

  std::vector<SATag> getSATag();
  std::vector<std::string> splitText(std::string text, char delimiter);
   std::vector<ReadParser::Cigar> getCigarByString(std::string cigartext);
   bool isProperlyAligned();
  
};

#endif