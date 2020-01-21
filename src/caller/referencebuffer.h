#ifndef REFERENCEBUFFER_H
#define REFERENCEBUFFER_H
#include "samplestat.h"
#include "filemanager.h"
#include <stdint.h>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <fasta/fastareader.h>

class ReferenceBuffer
{
private:
FastaReader fastareader;
  struct SequenceData
  {
    int32_t MaxPos = 0;
    int32_t MaxEnd = 0;
    std::string Sequence;
  };
  std::vector<SequenceData> datalist;
  std::string getSequenceData(std::string chr, int32_t pos, int32_t end);
  bool isCorrectRange(int32_t pos, int32_t end);

public:
  ReferenceBuffer();
  void setFastaReader(FastaReader fasta);
  void addReferenceBuffer(std::string chr, int32_t pos, int32_t end);
};

#endif