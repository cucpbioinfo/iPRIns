#ifndef FASTAREADER_H
#define FASTAREADER_H
#include <string>
#include <vector>
#include <map>
#include <stdint.h>
#include "faidx.h"
#include "fasta.h"
#include <algorithm>

class FastaReader
{
private:
  Faidx faidx;
  Fasta fasta;

public:
  // Fasta(std::string t_filepath, std::string t_indexfilepath);
  FastaReader();
  void setFilePath(std::string filepath);
  std::string getFilePath();
  bool hasFilePath();
  void setIndexFilePath(std::string indexfilepath);
  std::string getIndexFilePath();
  bool hasIndexFilePath();
  void exitIfNoFilePath();
  bool FileExists(std::string filename);
  void initialize();
  std::string getSeqbyPosition(std::string chromosome, int64_t start, int64_t end);

  void replaceSeqToUppercase(std::string *str)
  {
    std::transform((*str).begin(), (*str).end(), (*str).begin(), ::toupper);
  }
};

#endif