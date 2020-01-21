#ifndef REFININGSV_H
#define REFININGSV_H
#include "samplestat.h"
#include "filemanager.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refiningsv.h"
#include <fasta/fastareader.h>
#include "readparser.h"
#include "evidence.h"
#include "readdepthhelper.h"

class RefiningSV
{
private:
protected:
  bam1_t *read = NULL;
  bam_hdr_t *bam_header = NULL;
  samFile *inFile = NULL;
  hts_idx_t *bam_index = NULL;

  bool resultFromStart = false;

  SampleStat *samplestat;
  FileManager *filepath;
  FastaReader fastareader;
  ReadParser readparser;
  Evidence evidence;

public:
  struct AlternativeSA
  {
    std::string chr;
    int32_t pos=0;
  };

  struct MatchRead
  {
    std::string poschr;
    std::string endchr;
    int NumberOfMatchRead = 0;
    std::vector<int> MatchLists;
    std::vector<uint8_t> MapQLists;
    std::string Sequence;
    int maxMatchSequence = 0;
    int32_t maxSC = 0;
    int32_t maxAlterSC = 0;
    bool alignWithSoftClipped = false;
    std::vector<AlternativeSA> AltSA;
  };

public:
  RefiningSV();
  ~RefiningSV();
  Evidence variantresult;

  int getMaxIntFromVector(std::vector<int> value);
  uint8_t getMaxUInt8FromVector(std::vector<uint8_t> value);

  void setFilePath(FileManager *filepath);
  void prepareBamReader();
  void setHtsIndex(hts_idx_t *index);

  void setFastaReader(FastaReader fasta);
  void setSampleStat(SampleStat *samplepath_T);
  void setBamHeader(bam_hdr_t *bam_header);

  void setEvidence(Evidence evidence);

  void replaceSeqToUppercase(std::string *seq);

  std::string convertRangeToString(std::string chr, uint32_t pos, uint32_t end);
  void setRead(bam1_t *read);

  bool confirmAreaBySCAtStart(std::string chr, uint32_t pos, uint32_t end);
  bool confirmAreaBySCAtEnd(std::string chr, uint32_t pos, uint32_t end);

  bool isMatchRef(std::string chr, int32_t pos, int32_t end, std::string secondchr, int32_t secondpos, int32_t secondend);

  virtual void execute() = 0;

  std::string getResultVCFFormat();
  Evidence getVariantResult();
  int getReadDepthAtPosition(const char *range, int32_t pos);

  void calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *listPosition);

  bool haveIndel(std::vector<ReadParser::Cigar> cigar);
  int32_t getDivider(int value, int top, int down, int minimum);
};

#endif