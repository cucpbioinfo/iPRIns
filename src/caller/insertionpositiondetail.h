#ifndef INSERTIONPOSITIONDETAIL_H
#define INSERTIONPOSITIONDETAIL_H
#include "samplestat.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refiningsv.h"

class InsertionPositionDetail
{
private:
    int32_t position = 0;
    std::vector<uint8_t> mappingqualitylist;
    int longmapping = 0;
    std::vector<std::string> seqlist;

public:
    InsertionPositionDetail();
    void setPosition(int32_t position);
    int32_t getPosition();
    void addMapQ(uint8_t mapq);
    void addLongMapping(int longmapping);
    int getLongMapping();
    int getFrequency();
    uint8_t getMaxMapQ();
    uint8_t getMinMapQ();
    void addSeqList(std::string seq);
    std::vector<std::string> getSeqList();
    void setSeqList(std::vector<std::string> seqlist);

    std::vector<uint8_t> getMapQList();

    bool operator<(const InsertionPositionDetail &rhs) const
    {
        return (position < rhs.position);
    }
};

#endif