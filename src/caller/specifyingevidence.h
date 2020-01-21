#ifndef SPECIFYINGEVIDENCE_H
#define SPECIFYINGEVIDENCE_H
#include "readparser.h"
#include "evidence.h"
#include "readdepthhelper.h"

class SpecifyingEvidence
{
protected:
    int maxRemainBufferWriteEvidenceFile = 10;
    void writeBufferEvidenceFile();
  bam1_t *read;
  ReadParser readparser;
    ReadDepthHelper *readdepthHelper;
public:
    int getMaxRemainBufferWriteEvidenceFile() const;

    void setMaxRemainBufferWriteEvidenceFile(int bufferWriteEvidenceFile);

protected:
    bam_hdr_t *bam_header;
  SampleStat *samplestat;
  std::string outputpath;

  std::string svtype;

  std::vector<Evidence> preCollectSV;
  std::vector<Evidence> finalEvidence;
  
public:
  SpecifyingEvidence();
  bool checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped);
  int findOverlapped(int32_t overlapped, int32_t pos, int32_t mpos);
  int findOverlappedOnlyPos(uint32_t overlapped, uint32_t pos);
  void setRead(bam1_t *read, bam_hdr_t *bam_header);
  void setSampleStat(SampleStat *m_samplestat);
  void showAllFinalEvidence();
  void showSizeFinalEvidence();
  virtual void updateRead() = 0;
  void setOutputPath(std::string path);
  void writeFinalEvidenceAndClear();
  
  

    void setReadDepthHelper(ReadDepthHelper *readdepthHelper);
};

#endif