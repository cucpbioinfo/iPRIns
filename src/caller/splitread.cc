#include "splitread.h"
#include <iostream>
#include <fstream>
SplitRead::SplitRead(std::string chrname, ReadParser *readparser, SampleStat *samplestat, FileManager *filepath)
{
    SplitRead::readparser = readparser;
    SplitRead::samplestat = samplestat;
    SplitRead::chrname = chrname;
    SplitRead::filepath = filepath;
}

void SplitRead::updateRead()
{
    findInsertionInRead();
    satag = readparser->getSATag();
    if (satag.size() == 0)
    {
        return;
    }

    findInsertion();
}

void SplitRead::findInsertionInRead()
{
    std::vector<ReadParser::Cigar> cigar = readparser->getCigar();

    if (cigar.size() <= 2)
    {
        return;
    }

    // CIGAR = 20S71M6D43M54D117M = 20 + 71(M) + 43(M) + 117(M)
    int32_t incrementPos = 0;

    for (int i = 0; i < cigar.size(); i++)
    {
        if (cigar.at(i).getOperatorName() == 'M')
        {
            incrementPos += cigar.at(i).getLength();
        }

        if (cigar.at(i).getOperatorName() == 'I')
        {

            if (cigar.at(i).getLength() > 10)
            {

                mapSmallINS[std::make_pair(readparser->getPos() + incrementPos, readparser->getPos() + incrementPos + cigar.at(i).getLength())].NumberOfMatchRead++;
                mapSmallINS[std::make_pair(readparser->getPos() + incrementPos, readparser->getPos() + incrementPos + cigar.at(i).getLength())].MatchLists.push_back(cigar.at(i).getLength());
                mapSmallINS[std::make_pair(readparser->getPos() + incrementPos, readparser->getPos() + incrementPos + cigar.at(i).getLength())].MapQLists.push_back(readparser->getMapQuality());

            }

            incrementPos += cigar.at(i).getLength();
        }
    }
}

void SplitRead::findInsertion()
{

    for (ReadParser::SATag sa : satag)
    {

        if (sa.cigar.size() != 2)
        {
            continue;
        }

        if (readparser->getChromosomeNameString() != sa.chrname)
        {
            continue;
        }

        if (readparser->isReverse() && sa.strand == "-")
        {
            goto findInsertion;
        }

        if (!readparser->isReverse() && sa.strand == "+")
        {
            goto findInsertion;
        }

    findInsertion:

        // if (readparser->hasFirstCigarSoftclipped() && sa.cigar.at(0).getOperatorName() == 'S')
        // {

        //     if (checkBetween(readparser->getPos(), sa.pos, getDivider(samplestat->getReadLength(), 1, 10, 1)))
        //     {
        //         mapINS[std::make_pair(readparser->getPos(), sa.pos)].NumberOfMatchRead++;
        //         mapINS[std::make_pair(readparser->getPos(), sa.pos)].MatchLists.push_back(sa.cigar.at(sa.cigar.size() - 1).getLength());
        //         mapINS[std::make_pair(readparser->getPos(), sa.pos)].MapQLists.push_back(readparser->getMapQuality());
        //     }
        // }

        if (readparser->hasLastCigarSoftclipped() && sa.cigar.at(0).getOperatorName() == 'S')
        {
            int softclip = readparser->getSoftClippedSequenceEnd().size();
            int match = 0;
            if (sa.cigar.at(sa.cigar.size() - 1).getOperatorName() == 'M')
            {
                match = sa.cigar.at(sa.cigar.size() - 1).getLength();
            }

            if (match == 0)
            {
                continue;
            }

            int lsize = softclip - match;

            if (lsize<50)
            {
                continue;
            }

            if (checkBetween(readparser->getEnd(), sa.pos, getDivider(samplestat->getReadLength(), 1, 10, 1)))
            {
                mapINS[std::make_pair(readparser->getEnd(), sa.pos)].NumberOfMatchRead++;
                mapINS[std::make_pair(readparser->getEnd(), sa.pos)].MatchLists.push_back(lsize);
                mapINS[std::make_pair(readparser->getEnd(), sa.pos)].MapQLists.push_back(readparser->getMapQuality());
            }
        }
    }
}


void SplitRead::printResult()
{
    printInsertion();
}

std::vector<Evidence> SplitRead::convertMapToEvidenceList(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *mapSV, std::string svtype, std::string mark)
{
    std::vector<Evidence> vecTemp;
    for (auto const &x : *mapSV)
    {
        if (x.second.MapQLists.size() >= 1)
        {
            Evidence evidence;
            evidence.setPos(x.first.first);
            evidence.setEnd(x.first.second);
            evidence.setChr(chrname);
            evidence.setEndChr(chrname);
            evidence.setFrequency(x.second.MapQLists.size());
            evidence.setVariantType(svtype);
            evidence.setMark(mark);
            evidence.setMapQList(x.second.MapQLists);

            vecTemp.push_back(evidence);
        }
    }

    return vecTemp;
}

void SplitRead::mergeEvidence(std::vector<Evidence> *vecTemp)
{
    std::vector<Evidence> newEvidenceTempList;
    std::sort(vecTemp->begin(), vecTemp->end());

    for (auto m : *vecTemp)
    {
        bool incremented = false;
        for (int32_t i = 0; i < newEvidenceTempList.size(); i++)
        {
            if (checkBetween(m.getPos(), newEvidenceTempList.at(i).getPos(), 2) && checkBetween(m.getEnd(), newEvidenceTempList.at(i).getEnd(), 2))
            {
                incremented = true;
                for (auto x : *m.getMapQVector())
                {
                    newEvidenceTempList.at(i).addMapQ(x);
                    newEvidenceTempList.at(i).incrementFrequency();
                }
                break;
            }
        }

        if (!incremented)
        {
            newEvidenceTempList.push_back(m);
        }
    }

    *vecTemp = newEvidenceTempList;
}

void SplitRead::setAllCIEvidence(std::vector<Evidence> *elist, int32_t rangePos)
{
    for (int32_t i = 0; i < elist->size(); i++)
    {
        elist->at(i).setCiPosLeft(-rangePos);
        elist->at(i).setCiPosRight(rangePos);
        elist->at(i).setCiEndLeft(-rangePos);
        elist->at(i).setCiEndRight(rangePos);
    }
}

void SplitRead::filterFrequencyLowerThan(int number, std::vector<Evidence> *elist)
{
    std::vector<Evidence> newEvidenceTempList;

    for (auto n : *elist)
    {
        if (n.getFrequency() <= number)
        {
            continue;
        }

        newEvidenceTempList.push_back(n);
    }

    *elist = newEvidenceTempList;
}

void SplitRead::filterMapQLowerThan(uint8_t mapq, std::vector<Evidence> *elist)
{
    std::vector<Evidence> newEvidenceTempList;

    for (auto n : *elist)
    {
        if (n.getMaxMapQ() < mapq)
        {
            continue;
        }

        newEvidenceTempList.push_back(n);
    }

    *elist = newEvidenceTempList;
}

void SplitRead::filterEvidenceList(std::vector<Evidence> *elist)
{
    std::vector<Evidence> newEvidenceTempList;

    for (auto n : *elist)
    {
        if (n.getMaxMapQ() < 60)
        {
            continue;
        }

        if (n.getMinMapQ() == 0)
        {
            continue;
        }

        newEvidenceTempList.push_back(n);
    }

    *elist = newEvidenceTempList;
}

void SplitRead::filterLengthMinEvidenceList(std::vector<Evidence> *elist, int32_t min)
{
    std::vector<Evidence> newEvidenceTempList;

    for (auto n : *elist)
    {
        if (n.getSvLength() < min)
        {
            continue;
        }

        newEvidenceTempList.push_back(n);
    }

    *elist = newEvidenceTempList;
}

void SplitRead::filterLengthMaxEvidenceList(std::vector<Evidence> *elist, int32_t max)
{
    std::vector<Evidence> newEvidenceTempList;

    for (auto n : *elist)
    {
        if (n.getSvLength() > max)
        {
            continue;
        }

        newEvidenceTempList.push_back(n);
    }

    *elist = newEvidenceTempList;
}


void SplitRead::printInsertion()
{
    auto vecTemp = convertMapToEvidenceList(&mapINS, "INS", "SR");
    // filterLengthMinEvidenceList(&vecTemp, 49);
    // setAllCIEvidence(&vecTemp, samplestat->getReadLength() * 2);
    filterFrequencyLowerThan(0, &vecTemp);
    for (auto x : vecTemp)
    {
        writeFile(x);
    }
}

void SplitRead::printSmallInsertion()
{
    auto vecTemp = convertMapToEvidenceList(&mapSmallINS, "INS", "SINS");
    // mergeEvidence(&vecTemp);

    filterLengthMinEvidenceList(&vecTemp, 49);
    filterFrequencyLowerThan(1, &vecTemp);

    for (auto x : vecTemp)
    {
        writeFile(x);
    }
}


std::vector<Evidence> SplitRead::convertMapToEvidenceListTRA(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *mapSV, std::string svtype, std::string mark)
{
    std::vector<Evidence> vecTemp;
    for (auto const &x : *mapSV)
    {
        if (x.second.MapQLists.size() >= 1)
        {
            Evidence evidence;
            evidence.setPos(x.first.first);
            evidence.setEnd(x.first.second);
            evidence.setChr(x.second.poschr);
            std::cout << x.second.poschr << std::endl;
            evidence.setEndChr(x.second.endchr);
            evidence.setFrequency(x.second.MapQLists.size());
            evidence.setVariantType(svtype);
            evidence.setMark(mark);
            evidence.setMapQList(x.second.MapQLists);

            vecTemp.push_back(evidence);
        }
    }

    return vecTemp;
}

void SplitRead::removeDuplicateResult(std::vector<Evidence> *vec)
{
    std::vector<Evidence> temp;
    bool added = false;
    for (auto n : *vec)
    {
        added = false;
        for (auto m : temp)
        {
            if (n.getPos() == m.getPos() && n.getEnd() == m.getEnd())
            {
                continue;
            }

            if (checkBetween(n.getPos(), m.getPos(), 200) && checkBetween(n.getEnd(), m.getEnd(), 200))
            {
                added = true;
                break;
            }
        }

        if (!added)
        {
            temp.push_back(n);
        }
    }

    *vec = temp;
}

bool SplitRead::checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped)
{

    if (targetPos - overlapped > pos)
    {
        return false;
    }

    if (targetPos + overlapped < pos)
    {
        return false;
    }

    return true;
}

int SplitRead::writeFile(Evidence vr)
{
    std::ofstream myfile;
    myfile.open(filepath->getOutputPath() + "/analysis/splitread/" + vr.getChr() + "." + vr.getVariantType() + ".txt", std::ios_base::app);
    vr.setID("iPRIns" + std::to_string(vcfIdNumber));
    myfile << vr.getResultVcfFormatString() << std::endl;
    vcfIdNumber++;
    myfile.close();
    return 0;
}

int SplitRead::getDivider(int value, int top, int down, int minimum)
{
    auto returnvalue = (int)(float(value) * (float(top) / float(down)));

    if (returnvalue > minimum)
    {
        return returnvalue;
    }

    return minimum;
}