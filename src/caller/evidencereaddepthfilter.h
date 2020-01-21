#ifndef EVIDENCEREADDEPTHFILTER_H
#define EVIDENCEREADDEPTHFILTER_H


#include "filemanager.h"
#include "../fasta/fastareader.h"

class EvidenceReadDepthFilter {
private:
    FileManager *filepath;
    FastaReader fastaReader;
public:
    EvidenceReadDepthFilter(FileManager *filepath);
    bool haveNoneSymbol(std::string chr,int32_t start, int32_t end);
};

#endif