#include "samplestat.h"
#include <iostream>
#include <cmath>
#include <iterator>
#include <map>

void SampleStat::setSamplePath(std::string t_sample_path)
{
    sample_path = t_sample_path;
}

SampleStat::SampleStat()
{
    countMax = 100000;
}

void SampleStat::setNumberOfRead(int number)
{
    countMax = number;
}

SampleStat::SampleStat(std::string samplepath)
{
    sample_path = samplepath;
}

void SampleStat::findReadLength()
{
    samFile *fp_in = hts_open(sample_path.c_str(), "r"); //open bam file
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);             //read header
    bam1_t *aln = bam_init1();                           //initialize an alignment

    int count = 0;

    std::map<int, int> mapReadLength;
    while (sam_read1(fp_in, bamHdr, aln) > 0)
    {
        mapReadLength[aln->core.l_qseq]++;
        count++;
        if (count > 10000)
        {
            break;
        }
    }

    int maxhit = 0;
    int32_t lengthSeq = 0;
    for (const auto &p : mapReadLength)
    {
        if (maxhit < p.second)
        {
            maxhit = p.second;
            lengthSeq = p.first;
        }
    }

    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(fp_in);

    read_length = lengthSeq;
}

int32_t SampleStat::findSDSampleStat(std::vector<int32_t> *insertlist, int32_t insertsize_avg, int32_t limitsize, int32_t minsize)
{
    int count = 0;
    int64_t UPPER = 0;

    for (auto n : *insertlist)
    {
        if (n > limitsize)
        {
            continue;
        }

        if (n < minsize)
        {
            continue;
        }

        UPPER += pow((n - (insertsize_avg)), 2);
        count++;
    }

    return (int32_t)sqrt(UPPER / (int64_t)(count));
}

int32_t SampleStat::findMedianSampleStat(std::vector<int32_t> *insertlist, int32_t limitsize, int32_t minsize)
{
    int64_t sumINS = 0;
    int64_t count = 0;
    for (auto n : *insertlist)
    {
        if (n > limitsize)
        {
            continue;
        }

        if (n < minsize)
        {
            continue;
        }

        sumINS += n;
        count++;
    }

    return (int32_t)(sumINS / count);
}

void SampleStat::execute()
{
    SampleStat();
    findReadLength();
    auto insertsizelist = getInsertSizeList(countMax);

    int32_t avginsert = findMedianSampleStat(&insertsizelist, 100000, 0);
    int32_t sdinsert = findSDSampleStat(&insertsizelist, avginsert, 100000, 0);

    avginsert = findMedianSampleStat(&insertsizelist, avginsert + (sdinsert * 3), avginsert - (sdinsert * 3));
    sdinsert = findSDSampleStat(&insertsizelist, avginsert, avginsert + (sdinsert * 3), avginsert - (sdinsert * 3));

    insertsize_avg = avginsert;
    insertsize_sd = sdinsert;
}

std::vector<int32_t> SampleStat::getInsertSizeList(int64_t numberofread)
{
    samFile *fp_in = hts_open(sample_path.c_str(), "r");
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    bam1_t *aln = bam_init1();

    int count = 0;
    int64_t sumINS = 0;
    std::string chr;
    std::vector<int32_t> insertsizeList;
    while (sam_read1(fp_in, bamHdr, aln) > 0)
    {
        chr = bamHdr->target_name[aln->core.tid];
        if (!(chr == "1" || chr == "chr1"))
        {
            continue;
        }

        if (count > countMax)
        {
            break;
        }

        if (!(aln->core.flag & BAM_FPROPER_PAIR))
        {
            continue;
        }

        if (aln->core.tid != aln->core.mtid)
        {
            continue;
        }

        if ((aln->core.flag & BAM_FREAD2))
        {
            continue;
        }

        if ((aln->core.flag & BAM_FREVERSE))
        {
            continue;
        }

        if (!(aln->core.flag & BAM_FMREVERSE))
        {
            continue;
        }

        int32_t pos = aln->core.pos + 1;
        int32_t matepos = aln->core.mpos + 1;

        if (pos == 0)
        {
            continue;
        }

        if (matepos == 0)
        {
            continue;
        }

        if (pos > matepos)
        {
            continue;
        }

        if (aln->core.qual < 30)
        {
            continue;
        }

        if ((matepos + read_length - pos) > 0 && (matepos + read_length - pos) < 100000)
        {
            // if (matepos - pos < 200)
            // {
            //     std::cout << matepos - pos << std::endl;
            // }

            insertsizeList.push_back(matepos + read_length - pos);
            count++;
            
        }
        else
        {
            continue;
        }
    }

    std::cout <<"count : " << count << std::endl;

    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(fp_in);

    return insertsizeList;
}

int32_t SampleStat::getReadLength()
{
    return read_length;
}

int32_t SampleStat::getAverageSampleStat()
{
    return insertsize_avg;
}

int32_t SampleStat::getSDSampleStat()
{
    return insertsize_sd;
}