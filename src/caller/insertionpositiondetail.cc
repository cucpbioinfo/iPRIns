#include "insertionpositiondetail.h"

InsertionPositionDetail::InsertionPositionDetail()
{ 
    
}

void InsertionPositionDetail::setPosition(int32_t position)
{
    InsertionPositionDetail::position = position;
}

int32_t InsertionPositionDetail::getPosition()
{
    return InsertionPositionDetail::position;
}

void InsertionPositionDetail::addMapQ(uint8_t mapq)
{
    mappingqualitylist.push_back(mapq);
}

void InsertionPositionDetail::addLongMapping(int longmapping)
{
    if (InsertionPositionDetail::longmapping < longmapping)
    {
        InsertionPositionDetail::longmapping = longmapping;
    }
}

std::vector<uint8_t> InsertionPositionDetail::getMapQList()
{
    return InsertionPositionDetail::mappingqualitylist;
}

int InsertionPositionDetail::getLongMapping()
{
    return InsertionPositionDetail::longmapping;
}

int InsertionPositionDetail::getFrequency()
{
    return mappingqualitylist.size();
}

void InsertionPositionDetail::addSeqList(std::string seq)
{
    seqlist.push_back(seq);    
}

std::vector<std::string> InsertionPositionDetail::getSeqList()
{
    // std::cout << "seqlist : " << seqlist.size() << std::endl;
    // std::cout << "seq : " << seqlist[0] << std::endl;
    return seqlist;
}

void InsertionPositionDetail::setSeqList(std::vector<std::string> seqlist)
{
    InsertionPositionDetail::seqlist = seqlist;
}

uint8_t InsertionPositionDetail::getMaxMapQ()
{
    uint8_t max = 0;
    for (auto n : mappingqualitylist)
    {
        if (n > max)
        {
            max = n;
        }
    }

    return max;
}

uint8_t InsertionPositionDetail::getMinMapQ()
{
    uint8_t min = 255;
    for (auto n : mappingqualitylist)
    {
        if (n < min)
        {
            min = n;
        }
    }

    if (min == 255)
    {
        return 0;
    }

    return min;
}