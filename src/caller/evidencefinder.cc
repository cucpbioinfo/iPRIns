#include "evidencefinder.h"
#include <stdio.h>
#include <string>
#include <iostream>
#include <cstring>
#include <vector>
#include <iterator>
#include <algorithm>
#include "specifyingevidenceinsertion.h"
#include <unistd.h>
#include "readdepthhelper.h"
#include "splitread.h"

EvidenceFinder::EvidenceFinder(SampleStat *samplestat_T, FileManager *filepath_T, std::string *target_chromosome_T)
{
    samplestat = samplestat_T;
    filepath = filepath_T;
    target_chromosome = target_chromosome_T;
}

void getRound(uint32_t *x, uint32_t *max, uint32_t *round)
{
    *round = (*x / *max) * *max;
}

void EvidenceFinder::setHtsIndex(hts_idx_t *index)
{
    bam_index = index;
}

void EvidenceFinder::findEvidence()
{
    hts_itr_t *iterT = NULL;
    samFile *inT = NULL;
    // bam_hdr_t *bam_header = NULL;
    inT = sam_open(filepath->getSamplePath().c_str(), "r");
    if (inT == NULL)
        return;

    // target_chromosome std::string to array char
    int n = target_chromosome->length();
    char char_array_chrTarget[n + 1];
    strcpy(char_array_chrTarget, target_chromosome->c_str());
    iterT = sam_itr_querys(bam_index, bam_header, char_array_chrTarget);

    read = bam_init1();
    setupReadParser(read);
    int countT = 0;

    ReadDepthDetail readdepthdetail;
    uint32_t coverage = 0;
    uint32_t configRound = 250;
    uint32_t roundedPos = 0;
    uint32_t roundedCurrentPos = 0;

    int svtype;


    SpecifyingEvidenceInsertion seInsertion;
    seInsertion.setSampleStat(samplestat);
    seInsertion.setRead(read, bam_header);
    seInsertion.setOutputPath(filepath->getTempEvidencePath() + "/" + *target_chromosome + ".INS.txt");

    std::vector<ReadParser::Cigar> cigar;
    SplitRead splitread(*target_chromosome, &readparser, samplestat, filepath);
    int32_t markduppos = 0;
    int32_t markdupend = 0;
    int32_t markdupmatepos = 0;
    int32_t markdupmateend = 0;
    uint8_t markdupmapq = 0;

    while (sam_itr_next(inT, iterT, read) >= 0)
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

        if (markduppos==readparser.getPos() && markdupend==readparser.getEnd()) {
            if (markdupmatepos==readparser.getMatePos() && markdupmateend==readparser.getMateEnd()) {
                if (markdupmapq == readparser.getMapQuality()) {
                    continue;
                }
            }

        }

        markduppos = readparser.getPos();
        markdupend = readparser.getEnd();
        markdupmatepos = readparser.getMatePos();
        markdupmateend = readparser.getMateEnd();
        markdupmapq = readparser.getMapQuality();


        splitread.updateRead();

        currentPos = read->core.pos + 1;
        currentMPos = read->core.mpos + 1;

        cigar = readparser.getCigar();

        seInsertion.updateRead();

        //Read depth
        getRound(&currentPos, &configRound, &roundedPos);
        if (roundedCurrentPos == roundedPos)
        {
            readdepthdetail.RD++;

            if (cigar.at(0).getOperatorName() == 'S' && cigar.at(0).getLength() >= 8)
            {
                readdepthdetail.SCF++;
            }

            if (cigar.at(cigar.size() - 1).getOperatorName() == 'S' && cigar.at(cigar.size() - 1).getLength() >= 8)
            {
                readdepthdetail.SCL++;
            }

        }
        else
        {
            coverage = readdepthdetail.RD * samplestat->getReadLength() / configRound;
            ReadDepthLineSegment[roundedCurrentPos].RD += coverage;

            ReadDepthLineSegment[roundedCurrentPos].SCF += readdepthdetail.SCF;
            ReadDepthLineSegment[roundedCurrentPos].SCL += readdepthdetail.SCL;

            readdepthdetail.RD = 0;
            coverage = 0;

            readdepthdetail.SCF = 0;
            readdepthdetail.SCL = 0;

            roundedCurrentPos = roundedPos;
        }
    }

    splitread.printResult();

    std::cout << "✓ " << *target_chromosome << std::endl;

    hts_itr_destroy(iterT);
    bam_destroy1(read);
    sam_close(inT);

    ReadDepthHelper rdh(filepath);
    rdh.setRange(250);
    rdh.setReadDepthMap(ReadDepthLineSegment);
    rdh.setTargetChromosome(*target_chromosome);
    rdh.writeReadDepthLineFile(filepath->getReadDepthPath() + "/" + *target_chromosome + ".txt");
    rdh.findEvidence();
    rdh.calculateAvgReaddepth();
    rdh.writeReadDepthStat(filepath->getReadDepthStatPath());

    seInsertion.setReadDepthHelper(&rdh);
    seInsertion.done();
}


void EvidenceFinder::execute()
{
    findEvidence();
}

void EvidenceFinder::setSampleStat(SampleStat *samplestat_T)
{
    samplestat = samplestat_T;
}

void EvidenceFinder::setFilePath(FileManager *filepath_T)
{
    filepath = filepath_T;
}

void EvidenceFinder::setTargetChromosome(std::string *target_chromosome_T)
{
    target_chromosome = target_chromosome_T;
}

void EvidenceFinder::setBamHeader(bam_hdr_t *t_bam_header)
{
    bam_header = t_bam_header;
}

void EvidenceFinder::setupReadParser(bam1_t *alnT)
{
    read = alnT;
    readparser.setBamRead(alnT);
    readparser.setBamHeader(bam_header);
}