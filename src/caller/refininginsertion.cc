#include "refininginsertion.h"
#include <stdlib.h>
#include "smithwaterman.h"

RefiningInsertion::RefiningInsertion()
{
    variantresult.setVariantType("INS");
}

void RefiningInsertion::execute()
{
    variantresult.setChr(evidence.getChr());
    variantresult.setEndChr(evidence.getEndChr());
    prepareBamReader();

    if (evidence.getMark() == "SINS")
    {
        variantresult = evidence;
        variantresult.setQuailtyPass(true);
        return;
    }

    if (evidence.getMark() == "SR")
    {
        variantresult = evidence;
        variantresult.setQuailtyPass(true);
        return;
    }

    first();
}

void RefiningInsertion::first()
{
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getPos()+evidence.getCiPosLeft(),
                                                 evidence.getEnd()+evidence.getCiEndRight());

//    if (evidence.getPos() + evidence.getCiPosRight() - evidence.getPos() - evidence.getCiPosLeft() < 150)
//    {
//        return;
//    }

//     std::cout << "findRange : " << findRange << "\n" << evidence.getPos() + evidence.getCiPosLeft() - (evidence.getPos() + evidence.getCiPosRight()) << std::endl;
//    std::cout << ">>> : " << " " <<  evidence.getPos() << "," << readparser.getLengthSequence() << "," << evidence.getEnd() << "," << readparser.getLengthSequence() << std::endl;

    const char *range = findRange.c_str();
    refineStartToEnd(range);

    RefiningInsertion::convertMapSC();
    RefiningInsertion::clearMapSC();
    RefiningInsertion::findBreakpoint();
    RefiningInsertion::filterBreakpoint();

    if (variantresult.isQuailtyPass() == false && evidence.getMark() != "MATEUNMAPPED")
    {
        evidence.setMark("UNMERGE");
        RefiningInsertion::findUnmergeBreakpoint();
        RefiningInsertion::filterBreakpoint();
    }
}

void RefiningInsertion::refinewithReference()
{
}

void RefiningInsertion::refineStartToEnd(const char *range)
{
    hts_itr_t *iter = NULL;

    iter = sam_itr_querys(bam_index, bam_header, range);
    if (iter == NULL)
        return;
    read = bam_init1();
    readparser.setBamHeader(bam_header);
    readparser.setBamRead(read);

    int32_t startDist = evidence.getPos() + evidence.getCiPosRight();
    int32_t endDist = evidence.getPos() + evidence.getCiPosLeft();

    while (sam_itr_next(inFile, iter, read) >= 0)
    {

        if (readparser.isUnmapped())
        {
            continue;
        }

        if (readparser.isNotPassingFilters())
        {
            continue;
        }

        if (readparser.isPCR())
        {
            continue;
        }

        if (readparser.isSupplementaryAlignment())
        {
            continue;
        }

        // 10000
        if (readparser.getPos() < startDist - 10000)
        {
            break;
        }

        // 20000
        if (readparser.getPos() > endDist + 10000)
        {
            break;
        }

        auto cigar = readparser.getCigar();
        if (cigar.size() <= 1)
        {
            continue;
        }

        if (cigar.at(cigar.size() - 1).getOperatorName() == 'S' && cigar.at(cigar.size() - 1).getLength() >= 4)
        {
            mapSCEnd[readparser.getEnd()].addMapQ(readparser.getMapQuality());
            mapSCEnd[readparser.getEnd()].addLongMapping(cigar.at(cigar.size() - 1).getLength());
            mapSCEnd[readparser.getEnd()].setPosition(readparser.getEnd());
            mapSCEnd[readparser.getEnd()].addSeqList(readparser.getEdgeSeqFromEndSeq(cigar.at(cigar.size() - 1).getLength()));
        }

        if (cigar.at(0).getOperatorName() == 'S' && cigar.at(0).getLength() >= 4)
        {
            mapSCStart[readparser.getPos()].addMapQ(readparser.getMapQuality());
            mapSCStart[readparser.getPos()].addLongMapping(cigar.at(0).getLength());
            mapSCStart[readparser.getPos()].setPosition(readparser.getPos());
            mapSCStart[readparser.getPos()].addSeqList(readparser.getEdgeSeqFromStartSeq(cigar.at(0).getLength()));
        }
    }

    if (mapSCEnd.size() > 100)
    {
        return;
    }

    if (mapSCStart.size() > 100)
    {
        return;
    }

    hts_itr_destroy(iter);

    return;
}

void RefiningInsertion::filterBreakpoint()
{
    // std::cout << "# filterBreakpoint : " << vectorBP.size() << std::endl;
    std::sort(vectorBP.begin(), vectorBP.end());

    int maxFreq = 0;
    int maxLongMatch = 0;
    int maxscore = 0;

    for (BreakpointPosition n : vectorBP)
    {

        if (evidence.getMark() == "UNMERGE" && n.frequency >= 3)
        {

            variantresult.setPos(n.pos);
            variantresult.setEnd(n.end);
            variantresult.setFrequency(n.frequency);
            variantresult.setRPMapQ(*evidence.getMapQVector());
            variantresult.LNGMATCH = n.longmatch;

            variantresult.setMapQList(n.mappingqualitylist);
            variantresult.setChr(evidence.getChr());
            variantresult.setEndChr(evidence.getEndChr());
            variantresult.setQuailtyPass(true);
            variantresult.setMark(evidence.getMark());

            continue;
        }

        int32_t averagePos = 0;
        if (n.pos > n.end)
        {
            averagePos = n.pos;
        }
        else
        {
            averagePos = (n.end + n.pos) / 2;
        }

        if (n.frequency <= 2)
        {
            continue;
        }

        if (evidence.getMark() == "")
        {

            if (evidence.getFrequency() <= 3)
            {
                continue;
            }

            if (getMaxUInt8FromVector(n.mappingqualitylist) < 10)
            {
                continue;
            }

            if (getMaxUInt8FromVector(n.mappingqualitylist) < 60)
            {
                if (n.longmatch <= getDivider(samplestat->getReadLength(), 15, 100, 15))
                {
                    continue;
                }
                // continue;
            }
            else
            {
                if (n.longmatch <= getDivider(samplestat->getReadLength(), 20, 100, 15))
                {
                    continue;
                }
                // continue;
            }
        }
        else if (evidence.getMark() == "MATEUNMAPPED")
        {
            if (evidence.getFrequency() <= 4)
            {
                continue;
            }

            if (n.frequency <= 4)
            {
                continue;
            }

            if (getMaxUInt8FromVector(n.mappingqualitylist) < 20)
            {
                continue;
            }

            if (evidence.getMaxMapQ() < 40)
            {
                continue;
            }

            if (n.longmatch <= getDivider(samplestat->getReadLength(), 30, 100, 15))
            {
                continue;
            }
        }

        if (n.longmatch < getDivider(samplestat->getReadLength(), 15, 100, 15))
        {
            continue;
        }

        int score = (n.frequency) * (2 * n.longmatch);

        if (score <= maxscore)
        {
            continue;
        }

        maxscore = score;

        variantresult.setPos(averagePos);
        variantresult.setEnd(averagePos);
        variantresult.setFrequency(n.frequency);
        variantresult.setRPMapQ(*evidence.getMapQVector());
        variantresult.LNGMATCH = n.longmatch;

        variantresult.setMapQList(n.mappingqualitylist);
        variantresult.setChr(evidence.getChr());
        variantresult.setEndChr(evidence.getEndChr());
        variantresult.setQuailtyPass(true);
        variantresult.setMark(evidence.getMark());
    }
}


void RefiningInsertion::findUnmergeBreakpoint()
{
    vectorBP.clear();
    BreakpointPosition tempBP;
    tempBP.pos = (evidence.getPosDiscordantRead()+evidence.getLastPosDiscordantRead())/2;
    tempBP.end = (evidence.getPosDiscordantRead()+evidence.getLastPosDiscordantRead())/2;
    tempBP.frequency = evidence.getMapQVector()->size();
    tempBP.mappingqualitylist = *evidence.getMapQVector();
    vectorBP.push_back(tempBP);

//    if (tempBP.frequency <= 4)
//    {
//        continue;
//    }
//
//    added = true;
//    vectorBP.push_back(tempBP);
//
//
//    added = true;
//    vectorBP.push_back(tempBP)
//
//    for (InsertionPositionDetail n : vectorSCStart)
//    {
//        if (n.getLongMapping() < 15)
//        {
//            continue;
//        }
//
//        bool added;
//
//        std::vector<CountRefineSeq> mergeStart;
//        mergeStart = mergeString(n, false);
//
//
//
//        for (InsertionPositionDetail m : vectorSCEnd)
//        {
//
//            if (checkBetween(n.getPosition(), m.getPosition(),-50, 50))
//            {
//                BreakpointPosition tempBP;
//                tempBP.pos = n.getPosition();
//                tempBP.end = m.getPosition();
//
//                tempBP.frequency = n.getFrequency() + m.getFrequency();
//                tempBP.longmapstart = n.getLongMapping();
//                tempBP.longmapend = m.getLongMapping();
//
//                if (tempBP.frequency <= 4)
//                {
//                    continue;
//                }
//
//                added = true;
//                vectorBP.push_back(tempBP);
//            }
//
//        }
//    }
}

void RefiningInsertion::findBreakpoint()
{

    for (InsertionPositionDetail n : vectorSCStart)
    {
        if (n.getLongMapping() < 15)
        {
            continue;
        }

        bool added;

        // std::cout << n.getSeqList().size() << std::endl;
        // std::cout << "text :" ;
        // for (auto text:n.getSeqList())
        // {
        //     std::cout << text << std::endl;
        // }
        std::vector<CountRefineSeq> mergeStart;
        // if ((evidence.getMark() != "MATEUNMAPPED"))
        // {
        mergeStart = mergeString(n, false);
        // }



        for (InsertionPositionDetail m : vectorSCEnd)
        {

            if (m.getLongMapping() < 15)
            {
                continue;
            }

            // if (m.getFrequency() <= 1)
            // {
            //     continue;
            // }

            int frequency = 0;
            std::vector<uint8_t> mapq;
            int32_t longmatch = 0;
            // if ((evidence.getMark() != "MATEUNMAPPED"))
            // {
            std::vector<CountRefineSeq> mergeEnd = mergeString(m, true);
            std::vector<uint8_t> tempmapq;
            std::string seq1;
            std::string seq2;
            bool passoverlapped = getOverlappedSeq(mergeStart, mergeEnd, &frequency, &longmatch, &tempmapq, &seq1, &seq2);
            mapq = tempmapq;

            if (!passoverlapped)
            {
                continue;
            }

            if (checkBetween(n.getPosition(), m.getPosition(), -samplestat->getReadLength(), samplestat->getReadLength()))
            {
                BreakpointPosition tempBP;
                tempBP.pos = n.getPosition();
                tempBP.end = m.getPosition();
                tempBP.frequency = frequency;
                tempBP.longmatch = longmatch;
                tempBP.seq1 = seq1;
                tempBP.seq2 = seq2;

                tempBP.score = n.getFrequency() + m.getFrequency();
                tempBP.longmapstart = n.getLongMapping();
                tempBP.longmapend = m.getLongMapping();
                tempBP.mappingqualitylist = mapq;

                if (tempBP.frequency <= 2)
                {
                    continue;
                }

                added = true;
                vectorBP.push_back(tempBP);
            }
        }
    }
}

bool RefiningInsertion::checkBetween(int32_t pos, int32_t targetPos, int32_t minusoverlapped, int32_t plusoverlapped)
{
    if (targetPos + minusoverlapped > pos)
    {
        return false;
    }

    if (targetPos + plusoverlapped < pos)
    {
        return false;
    }

    return true;
}

void RefiningInsertion::convertMapSC()
{
    vectorSCStart = convertMapSCToVector(mapSCStart);
    std::sort(vectorSCStart.begin(), vectorSCStart.end());

    vectorSCEnd = convertMapSCToVector(mapSCEnd);
    std::sort(vectorSCEnd.begin(), vectorSCEnd.end());
}

std::vector<InsertionPositionDetail> RefiningInsertion::convertMapSCToVector(std::map<int32_t, InsertionPositionDetail> mapSC)
{

    std::vector<InsertionPositionDetail> vectorSC;
    int count = 0;
    for (auto x : mapSC)
    {
        vectorSC.push_back(x.second);
        // vectorSC->at(count).setSeqList(x.second.getSeqList());
        // std::cout << x.second.getSeqList()[0] << std::endl;
        // std::cout <<  vectorSC.at(count).getSeqList()[0] << std::endl;

        count++;
    }

    return vectorSC;
}

void RefiningInsertion::clearMapSC()
{
    mapSCStart.clear();
    mapSCEnd.clear();
}

std::vector<RefiningInsertion::CountRefineSeq> RefiningInsertion::mergeString(InsertionPositionDetail fragmentlist, bool fromstart)
{
    // compareEditDistance("ACCCCCACAGCTGTTACCCAGCGCCACACACAGAGCAGACGCTGAATCACTGCTTATTGACTGAATCAGCA", "ACCCCCACAGCTGTTACCCAGCGCCACACACAGAGCAGACGCTGAATCACTGCTTATTGACTGAATCAGCAATGGGGTACCT", false);

    // compareEditDistance("AATCACTGCTTATTGACTGAATCAGCAATGGGGT", "GCCACACACAGAGCAGACGCTGAATCACTGCTTATTGACTGAATCAGCAATGGGGT", false);
    // compareEditDistance("GCCACACACAGAGCAGACGCTGAATCACTGCTTATTGACTGAATCAGCAATGGGGT", "AATCACTGCTTATTGACTGAATCAGCAATGGGGT", false);
    // compareEditDistance("GCCACACACAGAGCAGACGCTGAATCACT", "GCCACACACAGAGCAGACGC", true);
    // compareEditDistance("GCCACACACAGAGCAGACGC", "GCCACACACAGAGCAGACGCTGAATCACT", true);

    std::vector<CountRefineSeq> tempSeq;

    // return tempSeq;

    int count = 0;
    for (std::string n : fragmentlist.getSeqList())
    //  for (int in=0;in < fragmentlist.)
    {

        bool added = false;
        for (int i = 0; i < tempSeq.size(); i++)
        {
            if (compareEditDistance(n, tempSeq.at(i).seq, fromstart))
            {
                if (n.size() > tempSeq.at(i).seq.size())
                {
                    tempSeq.at(i).seq = n;
                    tempSeq.at(i).count++;
                    tempSeq.at(i).mapqlist.push_back(fragmentlist.getMapQList().at(count));
                }
                added = true;
                // break;
            }
        }

        if (!added)
        {
            CountRefineSeq tempCRS;
            tempCRS.seq = n;
            tempCRS.count++;
            tempCRS.mapqlist.push_back(fragmentlist.getMapQList().at(count));
            tempSeq.push_back(tempCRS);
        }
        count++;
    }

    return tempSeq;
}

bool RefiningInsertion::compareEditDistance(std::string s1, std::string s2, bool fromstart)
{
    std::string temps1 = s1;
    std::string temps2 = s2;
    substringSeq(&temps1, &temps2, fromstart);

    EditDistance editdistance;
    int editpoint = editdistance.Compare(&temps1, &temps2);

    if (editpoint < 3)
    {
        return true;
    }

    return false;
}

void RefiningInsertion::substringSeq(std::string *s1, std::string *s2, bool fromstart)
{
    std::string temps1;
    std::string temps2;
    int diffsize = s1->size() - s2->size();
    if (s1->size() > s2->size())
    {
        if (!fromstart)
        {
            temps1 = s1->substr(s1->length() - s2->size());
            temps2 = *s2;
        }
        else
        {
            temps1 = s1->substr(0, s2->size());
            temps2 = *s2;
        }
    }
    else
    {
        if (!fromstart)
        {
            temps2 = s2->substr(s2->length() - s1->size());
            temps1 = *s1;
        }
        else
        {
            temps2 = s2->substr(0, s1->size());
            temps1 = *s1;
        }
    }

    *s1 = temps1;
    *s2 = temps2;
}

bool RefiningInsertion::getOverlappedSeq(std::vector<CountRefineSeq> startSeq, std::vector<CountRefineSeq> endSeq, int *frequency, int *longmatch, std::vector<uint8_t> *mapq, std::string *seq1, std::string *seq2)
{
    for (CountRefineSeq n : startSeq)
    {
        if (n.seq.size() < 20)
        {
            continue;
        }

        // std::cout << n.seq << " = " << n.count << " > " << n.seq.size() << std::endl;

        for (CountRefineSeq m : endSeq)
        {
            if (m.seq.size() < 20)
            {
                continue;
            }

            if (n.count + m.count <= 2)
            {
                continue;
            }

            if ((evidence.getMark() == ""))
            {
                SmithWaterman swm(&n.seq, 0, false);
                int maxmatch = swm.findMaxMatchInsertion(&m.seq);

                int seq1MatchSize = n.seq.size() - maxmatch;
                int seq2MatchSize = m.seq.size() - maxmatch;

                if (seq1MatchSize + seq2MatchSize + maxmatch < 50)
                {
                    continue;
                }

                if (maxmatch <= getDivider(samplestat->getReadLength(), 15, 100, 15))
                {
                    continue;
                }

                *frequency = n.count + m.count;
                std::vector<uint8_t> tempmapq;
                tempmapq.insert(tempmapq.end(), n.mapqlist.begin(), n.mapqlist.end());
                tempmapq.insert(tempmapq.end(), m.mapqlist.begin(), m.mapqlist.end());
                *mapq = tempmapq;
                *longmatch = maxmatch;

                if (maxmatch > 15)
                {
                    return true;
                }
            }
            else
            {

                if (n.seq.size() <= getDivider(samplestat->getReadLength(), 20, 100, 25))
                {
                    continue;
                }

                if (m.seq.size() <= getDivider(samplestat->getReadLength(), 20, 100, 25))
                {
                    continue;
                }

                *frequency = n.count + m.count;
                std::vector<uint8_t> tempmapq;
                tempmapq.insert(tempmapq.end(), n.mapqlist.begin(), n.mapqlist.end());
                tempmapq.insert(tempmapq.end(), m.mapqlist.begin(), m.mapqlist.end());
                *mapq = tempmapq;
                if (m.seq.size() >= n.seq.size())
                {
                    *longmatch = m.seq.size();
                }
                else
                {
                    *longmatch = n.seq.size();
                }

                return true;
            }
        }
    }

    return false;
}