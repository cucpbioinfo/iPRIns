#ifndef FAIDX_H
#define FAIDX_H
#include <string>
#include <vector>
#include <map>
#include <stdint.h>

class Faidx
{
private:
  std::string indexfilepath;

  struct FileFormat
  {
    std::string NAME;
    int64_t LENGTH;
    int64_t OFFSET;
    int64_t LINEBASES;
    int64_t LINEWIDTH;
    int64_t QUALOFFSET;
  };

  struct OffsetFile {
    int64_t start;
    int64_t end;
  };

  std::vector<FileFormat> faidxlist;

  int64_t getApproximateLineStart(std::string chromosome,int64_t start);
  int64_t getApproximateLineEnd(std::string chromosome,int64_t end);

  int64_t getApproximateOffsetStart(std::string chromosome,int64_t start);
  int64_t getApproximateOffsetEnd(std::string chromosome,int64_t end);
public:
  Faidx();
  void initialize();
  void setIndexFilePath(std::string indexfilepath);
  std::string getIndexFilePath();

  bool hasChromosomeName(std::string chr);
  int64_t getLengthbyChromosome(std::string chr);
  int64_t getOffsetbyChromosome(std::string chr);
  int64_t getLinebasesbyChromosome(std::string chr);
  int64_t getLinewidthbyChromosome(std::string chr);
  int64_t getQualoffsetbyChromosome(std::string chr);

  int64_t getOffsetEndByPosition(std::string chromosome, int64_t end);
  int64_t getOffsetStartByPosition(std::string chromosome, int64_t start);
};

#endif