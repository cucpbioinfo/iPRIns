#include "evidencereaddepthfilter.h"

EvidenceReadDepthFilter::EvidenceReadDepthFilter(FileManager *filepath)
{
    fastaReader.setFilePath((filepath->getReferencePath()));
    fastaReader.setIndexFilePath((filepath->getReferencePath()) + ".fai");
    fastaReader.initialize();
}

bool EvidenceReadDepthFilter::haveNoneSymbol(std::string chr,int32_t start, int32_t end)
{
    std::string seq  = fastaReader.getSeqbyPosition(chr,start,end);
    if (seq.find('N') != std::string::npos)
    {
        return true;
    }

    return false;
}