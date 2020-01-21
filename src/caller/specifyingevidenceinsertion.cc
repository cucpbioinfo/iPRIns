#include "specifyingevidenceinsertion.h"

SpecifyingEvidenceInsertion::SpecifyingEvidenceInsertion()
{
    svtype = "INS";
}

void SpecifyingEvidenceInsertion::updateRead()
{

    if (readparser.isUnmapped())
    {
        return;
    }

    currentPos = read->core.pos + 1;
    currentMPos = read->core.mpos + 1;

    if (readparser.isFirstRead() && !readparser.isReverse())
    {
        if (readparser.isMateUnmapped())
        {
            checkRange();
            return;
        }
    }

    if (readparser.isSecondRead() && readparser.isReverse())
    {
        if (readparser.isMateUnmapped())
        {
            checkRange();
            return;
        }
    }

    if (!readparser.isPairOnSameChromosome())
    {
        return;
    }

    if (readparser.isFirstRead())
    {
        if (readparser.isReverse())
        {
            return;
        }

        if (!readparser.isMateReverse())
        {
            return;
        }

        int diff = (readparser.getMatePos() + readparser.getLengthSequence()) - readparser.getPos();
        if (diff < 0)
        {
            return;
        }

        if (!(read->core.flag & BAM_FPROPER_PAIR))
        {
            return;
        }

        // if (diff < samplestat->getAverageSampleStat() - int32_t(1.654 * double(samplestat->getSDSampleStat())))
        if (diff < samplestat->getAverageSampleStat() - (2 * samplestat->getSDSampleStat()))
        {
            checkRange();
        }
        return;
    }
}

bool SpecifyingEvidenceInsertion::incrementSVFreq(int32_t overlappedpos, int32_t pos, int32_t mpos)
{
    bool added = false;
    for (int positionOverlapped = 0; positionOverlapped < preCollectSV.size(); positionOverlapped++)
    {
        if (readparser.isMateUnmapped())
        {
            if (readparser.isFirstRead())
            {
                if (checkBetween(pos, preCollectSV.at(positionOverlapped).getPosDiscordantRead(), overlappedpos))
                {
                    if (preCollectSV.at(positionOverlapped).getLastPosDiscordantRead() < currentPos)
                    {
                        preCollectSV.at(positionOverlapped).setLastPosDiscordantRead(currentPos);
                    }

                    preCollectSV.at(positionOverlapped).addMapQ(readparser.getMapQuality());

                    preCollectSV.at(positionOverlapped).incrementFrequency();
                    added = true;
                }
            }
            else
            {
                if (checkBetween(pos, preCollectSV.at(positionOverlapped).getPosDiscordantRead() + samplestat->getAverageSampleStat() + (3 * samplestat->getSDSampleStat()), overlappedpos))
                {
                    if (preCollectSV.at(positionOverlapped).getEnd() == 0)
                    {
                        preCollectSV.at(positionOverlapped).setEnd(currentPos);
                        preCollectSV.at(positionOverlapped).setEndDiscordantRead(currentPos);
                    }

                    if (preCollectSV.at(positionOverlapped).getLastEndDiscordantRead() < currentPos)
                    {
                        preCollectSV.at(positionOverlapped).setLastEndDiscordantRead(currentPos);
                    }

                    preCollectSV.at(positionOverlapped).addMapQ(readparser.getMapQuality());

                    preCollectSV.at(positionOverlapped).incrementFrequency();
                    added = true;
                }
            }
        }
        else
        {
            if (readparser.isSecondRead())
            {
                continue;
            }

            if (checkBetween(pos, preCollectSV.at(positionOverlapped).getPosDiscordantRead(), overlappedpos))
            {
                if (preCollectSV.at(positionOverlapped).getLastPosDiscordantRead() < currentPos)
                {
                    preCollectSV.at(positionOverlapped).setLastPosDiscordantRead(currentPos);
                }

                if (currentMPos != 0)
                {
                    if (preCollectSV.at(positionOverlapped).getEndDiscordantRead() > currentMPos)
                    {
                        preCollectSV.at(positionOverlapped).setEndDiscordantRead(currentMPos);
                    }

                    if (preCollectSV.at(positionOverlapped).getLastEndDiscordantRead() < currentMPos)
                    {
                        preCollectSV.at(positionOverlapped).setLastEndDiscordantRead(currentMPos);
                    }
                }

                preCollectSV.at(positionOverlapped).addMapQ(readparser.getMapQuality());

                preCollectSV.at(positionOverlapped).incrementFrequency();
                added = true;
            }
        }
    }

    return added;
}

void SpecifyingEvidenceInsertion::checkRange()
{
    bool added = false;
    int32_t merge = int32_t(samplestat->getAverageSampleStat()) + int32_t(samplestat->getSDSampleStat() * 3);
    // if (readparser.isMateUnmapped()) && readparser.isSecondRead())
    // {
    //     int32_t merge = int32_t(samplestat->getAverageSampleStat()) + int32_t(samplestat->getSDSampleStat() * 3);
    // }
    added = incrementSVFreq(merge, currentPos, currentMPos);

    checkProveEvidence();

    if (!added)
    {
        if (readparser.isSecondRead())
        {
            return;
        }

        Evidence evidence;
        evidence.setVariantType(svtype);
        evidence.setChr(readparser.getChromosomeNameString());
        evidence.setEndChr(readparser.getChromosomeNameString());

        if (readparser.isMateUnmapped())
        {
            evidence.setMark("MATEUNMAPPED");
            evidence.setPosDiscordantRead(readparser.getPos());
            evidence.setLastPosDiscordantRead(readparser.getPos());
        }
        else
        {
            evidence.setPosDiscordantRead(readparser.getPos());
            evidence.setLastPosDiscordantRead(readparser.getPos());
            evidence.setEndDiscordantRead(readparser.getMatePos());
            evidence.setLastEndDiscordantRead(readparser.getMatePos());
        }

        evidence.addMapQ(readparser.getMapQuality());
        evidence.incrementFrequency();

        preCollectSV.push_back(evidence);
    }
}

void SpecifyingEvidenceInsertion::proveEvidence(int index)
{
    int32_t plus = samplestat->getAverageSampleStat() + (3*samplestat->getSDSampleStat()) + samplestat->getReadLength();
    if (currentPos - plus > preCollectSV.at(index).getPosDiscordantRead())
    {
        if (filterEvidence(&preCollectSV.at(index)))
        {
            calculateVCF(&preCollectSV.at(index));
            finalEvidence.push_back(preCollectSV.at(index));
            preCollectSV.erase(preCollectSV.begin() + index);
            writeBufferEvidenceFile();
            preCollectSV.clear();
        }
        else
        {
            preCollectSV.erase(preCollectSV.begin() + index);
        }
    }
}

void SpecifyingEvidenceInsertion::calculateVCF(Evidence *evidence)
{

    int32_t pos = evidence->getPosDiscordantRead();

    int32_t lastpos = evidence->getLastPosDiscordantRead();
    int32_t end = evidence->getEndDiscordantRead();
    int32_t lastend = evidence->getLastEndDiscordantRead();

    if (evidence->getMark() == "MATEUNMAPPED")
    {
        // end = pos + samplestat->getAverageSampleStat() + samplestat->getSDSampleStat();
        evidence->setPos(pos);
        evidence->setEnd(end);

        evidence->setCiPosLeft(-samplestat->getReadLength());
        evidence->setCiPosRight(lastend-pos+samplestat->getReadLength());

        evidence->setCiEndLeft(pos - lastend);
        evidence->setCiEndRight(samplestat->getReadLength());
    }
    else
    {
        // end = pos + samplestat->getAverageSampleStat() + samplestat->getSDSampleStat();
        evidence->setPos(pos);
        evidence->setEnd(end);

        evidence->setCiPosLeft(-(2*samplestat->getReadLength()));
        evidence->setCiPosRight((lastend - pos)+(2*samplestat->getReadLength()));

        evidence->setCiEndLeft(pos - end-(2*samplestat->getReadLength()));
        evidence->setCiEndRight((2*samplestat->getReadLength()));
    }
    // else
    // {
    //     // end = pos + samplestat->getAverageSampleStat() + samplestat->getSDSampleStat();
    //     evidence->setPos(pos);
    //     evidence->setEnd(end);

    //     evidence->setCiPosLeft(-(samplestat->getReadLength()));
    //     evidence->setCiPosRight((lastend - pos)+(samplestat->getReadLength()));

    //     evidence->setCiEndLeft(pos - end-(samplestat->getReadLength()));
    //     evidence->setCiEndRight((samplestat->getReadLength()));
    // }
}

bool SpecifyingEvidenceInsertion::filterEvidence(Evidence *evidence)
{

    if (evidence->getFrequency() <= 2)
    {
        return false;
    }

    // if (evidence->getMaxMapQ() < 20)
    // {
    //     return false;
    // }

    if (evidence->getMark() == "MATEUNMAPPED")
    {
        if (evidence->getLastPosDiscordantRead()-evidence->getPosDiscordantRead()==0)
        {
            return false;
        }

        if (evidence->getLastEndDiscordantRead()-evidence->getEndDiscordantRead()==0)
        {
            return false;
        }

        if (evidence->getEnd() == 0)
        {
            return false;
        }
    }

    return true;
}

void SpecifyingEvidenceInsertion::checkProveEvidence()
{
    for (int i = 0; i < preCollectSV.size(); i++)
    {
        proveEvidence(i);
    }
}

void SpecifyingEvidenceInsertion::done()
{
    std::vector<Evidence> temp;
    for (auto n : finalEvidence)
    {
        if (readdepthHelper->isRangeDisorderByMorethanRD(n.getPosDiscordantRead(), n.getEndDiscordantRead(), readdepthHelper->getAvgReadDepth() + 30))
        {
            continue;
        }
        temp.push_back(n);
    }
    finalEvidence = temp;

    checkProveEvidence();
    writeFinalEvidenceAndClear();
}