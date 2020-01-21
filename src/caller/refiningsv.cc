#include "refiningsv.h"

RefiningSV::RefiningSV()
{
}



RefiningSV::~RefiningSV()
{

    if (bam_header != NULL)
    {
        bam_hdr_destroy(bam_header);
    }

    if (inFile != NULL)
    {
        sam_close(inFile);
    }

    if (read != NULL)
    {
        bam_destroy1(read);
    }

    // if (bam_index != NULL)
    // {
    //     hts_idx_destroy(bam_index);
    // }
}

int32_t RefiningSV::getDivider(int value, int top, int down, int minimum)
{
    auto returnvalue = (int32_t)(float(value) * (float(top) / float(down)));

    if (returnvalue > minimum)
    {
        return returnvalue;
    }

    return minimum;
}

void RefiningSV::setSampleStat(SampleStat *samplestat_T)
{
    samplestat = samplestat_T;
}

int RefiningSV::getReadDepthAtPosition(const char *range, int32_t pos)
{
    hts_itr_t *iter = NULL;

    iter = sam_itr_querys(bam_index, bam_header, range);
    if (iter == NULL)
        return false;
    read = bam_init1();
    readparser.setBamHeader(bam_header);
    readparser.setBamRead(read);

    std::map<int32_t, int> mapSCFirst;
    std::map<int32_t, int> mapSCLast;

    int count = 0;
    while (sam_itr_next(inFile, iter, read) >= 0)
    {

        if (readparser.isUnmapped())
        {
            continue;
        }

        if (readparser.getPos() <= pos && pos <= readparser.getPos() + samplestat->getReadLength())
        {
            count++;
        }
    }

    hts_itr_destroy(iter);

    return count;
}

int RefiningSV::getMaxIntFromVector(std::vector<int> value)
{
    int max = 0;
    for (auto n : value)
    {
        if (n > max)
        {
            max = n;
        }
    }

    return max;
}

uint8_t RefiningSV::getMaxUInt8FromVector(std::vector<uint8_t> value)
{
    uint8_t max = 0;
    for (auto n : value)
    {
        if (n > max)
        {
            max = n;
        }
    }

    return max;
}

void RefiningSV::setFastaReader(FastaReader fasta)
{
    fastareader = fasta;
}

void RefiningSV::setBamHeader(bam_hdr_t *t_bam_header)
{
    bam_header = t_bam_header;
    readparser.setBamHeader(bam_header);
}

void RefiningSV::setRead(bam1_t *alnT)
{
    read = alnT;
    readparser.setBamRead(alnT);
}

void RefiningSV::setEvidence(Evidence evidence_m)
{
    evidence = evidence_m;
}

void RefiningSV::setFilePath(FileManager *filepath_T)
{
    filepath = filepath_T;
}

std::string RefiningSV::convertRangeToString(std::string chr, uint32_t m_pos, uint32_t m_end)
{
    std::string pos = std::to_string(m_pos);
    std::string end = std::to_string(m_end);
    return chr + ":" + pos + "-" + end;
}

void RefiningSV::setHtsIndex(hts_idx_t *index)
{
    bam_index = index;
}

void RefiningSV::replaceSeqToUppercase(std::string *str)
{
    std::transform((*str).begin(), (*str).end(), (*str).begin(), ::toupper);
}

void RefiningSV::prepareBamReader()
{
    inFile = sam_open(filepath->getSamplePath().c_str(), "r");
    if (inFile == NULL)
        return;
    if ((bam_header = sam_hdr_read(inFile)) == 0)
        return;
    // bam_index = sam_index_load(inFile, filepath->getSamplePath()->c_str());
    // if (bam_index == NULL)
    //     return;
}

bool RefiningSV::confirmAreaBySCAtStart(std::string chr, uint32_t pos, uint32_t end)
{
    std::string findRange = convertRangeToString(chr, pos - samplestat->getReadLength(), end + samplestat->getReadLength());
    const char *range = findRange.c_str();
    hts_itr_t *iter = NULL;

    iter = sam_itr_querys(bam_index, bam_header, range);
    if (iter == NULL)
        return false;
    read = bam_init1();
    ReadParser readparser;
    readparser.setBamHeader(bam_header);
    readparser.setBamRead(read);

    int count = 0;

    while (sam_itr_next(inFile, iter, read) >= 0)
    {
        if (readparser.isUnmapped())
        {
            continue;
        }

        if (readparser.hasFirstCigarSoftclipped())
        {

            auto cigar = readparser.getCigar();

            uint32_t posSC = readparser.getPos() - cigar[0].getLength();
            uint32_t endSC = readparser.getPos();

            if (end <= endSC && pos >= posSC)
            {
                count++;
                // std::cout << readparser.getPosDiscordantRead() << " / "
                //           << readparser.getEndDiscordantRead() << " cigar : "
                //           << cigar[cigar.size() - 1].getLength() << " end real : "
                //           << readparser.getEndDiscordantRead() - cigar[cigar.size() - 1].getLength()
                //           << " "
                //           << readparser.getSequence()
                //           << std::endl;
            }
        }
    }

    hts_itr_destroy(iter);

    if (count >= 3)
    {
        return true;
    }

    return false;
}

bool RefiningSV::confirmAreaBySCAtEnd(std::string chr, uint32_t pos, uint32_t end)
{
    std::string findRange = convertRangeToString(chr, pos - samplestat->getReadLength(), end + samplestat->getReadLength());
    const char *range = findRange.c_str();
    hts_itr_t *iter = NULL;

    iter = sam_itr_querys(bam_index, bam_header, range);
    if (iter == NULL)
        return false;
    read = bam_init1();
    ReadParser readparser;
    readparser.setBamHeader(bam_header);
    readparser.setBamRead(read);

    int count = 0;

    while (sam_itr_next(inFile, iter, read) >= 0)
    {
        if (readparser.isUnmapped())
        {
            continue;
        }

        if (readparser.hasLastCigarSoftclipped())
        {

            auto cigar = readparser.getCigar();

            //Position before last softclipped
            uint32_t posSC = readparser.getEnd() - cigar[cigar.size() - 1].getLength();
            //Position after last softclipped
            uint32_t endSC = readparser.getEnd();

            if (posSC >= pos && posSC <= endSC)
            {
                //     count++;
                // std::cout << readparser.getPosDiscordantRead() << " / "
                //           << readparser.getEndDiscordantRead() << " cigar : "
                //           << cigar[cigar.size() - 1].getLength()
                //           << " end real : " << readparser.getEndDiscordantRead() - cigar[cigar.size() - 1].getLength()
                //           << " "
                //           << readparser.getSequence()
                //           << std::endl;
            }
        }
    }

    hts_itr_destroy(iter);

    if (count >= 3)
    {
        return true;
    }

    return false;
}

Evidence RefiningSV::getVariantResult()
{
    return variantresult;
}

std::string RefiningSV::getResultVCFFormat()
{
    return variantresult.getResultVcfFormatString();
}

// if compare 449280 = 449281 @ NGB
bool RefiningSV::isMatchRef(std::string chr, int32_t pos, int32_t end, std::string secondchr, int32_t secondpos, int32_t secondend)
{
    std::string referenceSeq = fastareader.getSeqbyPosition(chr, pos, end);
    // std::cout << "referenceSeq :" << referenceSeq << std::endl;

    std::string referenceSeqSecond = fastareader.getSeqbyPosition(secondchr, secondpos, secondend);
    // std::cout << "referenceSeqSecond :" << referenceSeqSecond << std::endl;

    replaceSeqToUppercase(&referenceSeq);
    replaceSeqToUppercase(&referenceSeqSecond);
    if (referenceSeq == referenceSeqSecond)
    {
        return true;
    }

    return false;
}

bool RefiningSV::haveIndel(std::vector<ReadParser::Cigar> cigar)
{
    for (auto n : cigar)
    {
        if (n.getOperatorName() == 'D')
        {
            return true;
        }

        if (n.getOperatorName() == 'I')
        {
            return true;
        }
    }

    return false;
}

void RefiningSV::calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *listPosition)
{

    int32_t bPos = 0;
    int32_t bEnd = 0;
    int32_t bHit = 0;
    uint8_t bMaxQuality = 0;
    int bMaxMatchSize = 0;
    int bFrequency = 0;
    std::vector<uint8_t> bMapQList;
    int32_t svlength = evidence.getEndDiscordantRead() - evidence.getPosDiscordantRead() - samplestat->getAverageSampleStat();
    // std::cout << "svlength :" << svlength << std::endl;
    int lastscore = 0;

    for (auto const &x : *listPosition)
    {
        int maxMatchSize = getMaxIntFromVector(x.second.MatchLists);

        uint8_t maxQuality = getMaxUInt8FromVector(x.second.MapQLists);

        if (maxMatchSize < 25)
        {
            continue;
        }

        if (maxMatchSize > 80)
        {
            continue;
        }

        if (x.second.maxAlterSC >= x.second.maxSC && x.second.alignWithSoftClipped)
        {
            continue;
        }

        bFrequency = x.second.NumberOfMatchRead;

        if (evidence.getSvLength() < 500)
        {
            if (bFrequency <= 2)
            {
                continue;
            }
            if (maxQuality == 0)
            {
                continue;
            }
        }
        else
        {

            if (bFrequency <= 1)
            {
                continue;
            }
        }

        if (isMatchRef(evidence.getChr(), x.first.first - 1, x.first.first + 20 - 1, evidence.getEndChr(), x.first.second, x.first.second + 20))
        {
            continue;
        }

        // confirm
        if (isMatchRef(evidence.getChr(), x.first.first - 1, x.first.first + 20 - 1, evidence.getEndChr(), x.first.second + 1, x.first.second + 20 + 1))
        {
            continue;
        }

        if (isMatchRef(evidence.getChr(), x.first.first - 1 + 2, x.first.first + 20 - 1 + 2, evidence.getEndChr(), x.first.second + 2, x.first.second + 20 + 2))
        {
            continue;
        }

        if (isMatchRef(evidence.getChr(), x.first.first + 1 - 20, x.first.first + 1, evidence.getEndChr(), x.first.second - 20, x.first.second))
        {
            continue;
        }

        if (isMatchRef(evidence.getChr(), x.first.first + 2 - 20, x.first.first + 1 + 2, evidence.getEndChr(), x.first.second + 20 + 2, x.first.second + 2))
        {
            continue;
        }

        int number = x.second.NumberOfMatchRead;

        int score = (number) * (2 * maxMatchSize);

        if (score > lastscore)
        {
            lastscore = score;
            bPos = x.first.first;
            bEnd = x.first.second;
            bHit = number;
            bMaxMatchSize = maxMatchSize;
            bMapQList = x.second.MapQLists;
        }
    }

    variantresult.setPos(bPos);
    variantresult.setEnd(bEnd);
    variantresult.setFrequency(bHit);
    variantresult.setMapQList(bMapQList);
    variantresult.setChr(evidence.getChr());
    variantresult.setEndChr(evidence.getEndChr());
    variantresult.LNGMATCH = bMaxMatchSize;

    if (bPos == 0)
    {
        return;
    }
    if (bEnd == 0)
    {
        return;
    }

    if (variantresult.getVariantType() == "DEL")
    {
        if (bPos == 0)
        {
            return;
        }
        if (bEnd == 0)
        {
            return;
        }

        if (bEnd - bPos > 50000)
        {
            return;
        }

        if (bEnd - bPos < 0)
        {
            return;
        }

        variantresult.setQuailtyPass(true);
        return;
    }

    if (variantresult.getVariantType() == "DUP")
    {
        variantresult.setQuailtyPass(true);
        return;
    }

    if (variantresult.getVariantType() == "INS")
    {
        variantresult.setQuailtyPass(true);
        return;
    }

    if (variantresult.getVariantType() == "INV")
    {
        variantresult.setQuailtyPass(true);
        return;
    }

    if (variantresult.getVariantType() == "BND")
    {
        variantresult.setQuailtyPass(true);
        return;
    }

    variantresult.setQuailtyPass(false);
}
