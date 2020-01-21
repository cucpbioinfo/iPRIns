#include "readparser.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>
#include <cstring>
#include <stdlib.h>
#include <ctype.h>
#include <sstream>
//#include <bits/stdc++.h>
void ReadParser::setBamRead(bam1_t *bam_read)
{
    source_bamread = bam_read;
}

void ReadParser::setBamHeader(bam_hdr_t *t_bam_header)
{
    bam_header = t_bam_header;
}

ReadParser::ReadParser()
{
}

ReadParser::~ReadParser()
{
}

std::vector<ReadParser::Cigar> ReadParser::getCigar()
{
    const auto cigar = bam_get_cigar(source_bamread);

    std::vector<ReadParser::Cigar> cigarlist;

    for (int k = 0; k < getLengthCigar(); k++)
    {

        const int op = bam_cigar_op(cigar[k]);
        const int ol = bam_cigar_oplen(cigar[k]);

        Cigar cigarL;
        cigarL.setOperatorName(op);
        cigarL.setLength(ol);
        cigarlist.push_back(cigarL);
    }

    return cigarlist;
}

int32_t ReadParser::getPos()
{
    return source_bamread->core.pos + 1;
}

int32_t ReadParser::getEnd()
{
    std::vector<ReadParser::Cigar> cigars = getCigar();
    int count = 0;
    int position = 0;
    int sizeCigar = cigars.size();
    for (Cigar c : cigars)
    {
        if (c.getOperatorName() == 'S' && position == 0)
        {
            count -= c.getLength();
        }
        else if (c.getOperatorName() == 'I')
        {
        }
        else if (c.getOperatorName() == 'H')
        {
        }
        else if (c.getOperatorName() == 'S' && position != sizeCigar)
        {
        }
        else
        {
            count += c.getLength();
        }
        position++;
    }

    return getPos() + count + 1;
}

int32_t ReadParser::getMatePos()
{
    return source_bamread->core.mpos + 1;
}

int32_t ReadParser::getMateEnd()
{
    return source_bamread->core.mpos + source_bamread->core.l_qseq + 1;
}

int32_t ReadParser::getInsertSize()
{
    if (getPos() < getMatePos())
    {
        return getMateEnd() - getPos();
    }

    return getEnd() - getMatePos();
}

int32_t ReadParser::getPosOfSeq()
{
    std::vector<ReadParser::Cigar> cigars = getCigar();

    int32_t shift = 0;
    if (cigars.at(0).getOperatorName() == 'S')
    {
        shift = (int32_t)cigars.at(0).getLength();
    }

    return getPos() - shift;
}

int32_t ReadParser::getEndOfSeq()
{
    return getPosOfSeq() + (int32_t)source_bamread->core.l_qseq;
}

bool ReadParser::isMateUnmapped()
{
    if ((source_bamread->core.flag & BAM_FMUNMAP))
    {
        return true;
    }

    return false;
}

bool ReadParser::isUnmapped()
{
    if ((source_bamread->core.flag & BAM_FUNMAP))
    {
        return true;
    }

    return false;
}

std::string ReadParser::getChromosomeNameString()
{
    std::string chr = bam_header->target_name[source_bamread->core.tid];
    return chr;
}

std::string ReadParser::getMateChromosomeNameString()
{
    std::string mchr = bam_header->target_name[source_bamread->core.mtid];
    return mchr;
}

bool ReadParser::isPairOnSameChromosome()
{
    if (bam_header->target_name[source_bamread->core.tid] == bam_header->target_name[source_bamread->core.mtid])
    {
        return true;
    }

    return false;
}

uint32_t ReadParser::getLengthCigar()
{
    return source_bamread->core.n_cigar;
}

void ReadParser::printCigar()
{
    const auto cigar = bam_get_cigar(source_bamread);

    for (int k = 0; k < getLengthCigar(); k++)
    {
        const int op = bam_cigar_op(cigar[k]);
        const int ol = bam_cigar_oplen(cigar[k]);

        // std::cout << op << "=" << ol << " ";
    }
}

bool ReadParser::hasLastCigarSoftclipped()
{
    const auto cigar = bam_get_cigar(source_bamread);
    int last = getLengthCigar() - 1;
    const int op = bam_cigar_op(cigar[last]);
    const int ol = bam_cigar_oplen(cigar[last]);

    if (op == BAM_CSOFT_CLIP)
    {
        return true;
    }

    return false;
}

bool ReadParser::hasFirstCigarSoftclipped()
{
    const auto cigar = bam_get_cigar(source_bamread);
    const int op = bam_cigar_op(cigar[0]);
    const int ol = bam_cigar_oplen(cigar[0]);

    if (op == BAM_CSOFT_CLIP)
    {
        return true;
    }

    return false;
}

bool ReadParser::hasLastCigarHardclipped()
{
    const auto cigar = bam_get_cigar(source_bamread);
    int last = getLengthCigar() - 1;
    const int op = bam_cigar_op(cigar[last]);
    const int ol = bam_cigar_oplen(cigar[last]);

    if (op == BAM_CHARD_CLIP)
    {
        return true;
    }

    return false;
}

bool ReadParser::hasFirstCigarHardclipped()
{
    const auto cigar = bam_get_cigar(source_bamread);
    const int op = bam_cigar_op(cigar[0]);
    const int ol = bam_cigar_oplen(cigar[0]);

    if (op == BAM_CHARD_CLIP)
    {
        return true;
    }

    return false;
}

std::string ReadParser::getSequence()
{
    std::string buf;
    for (int i = 0; i < source_bamread->core.l_qseq; ++i)
    {
        buf += seq_nt16_str[bam_seqi(bam_get_seq(source_bamread), i)];
    }
    return buf;
}

int *ReadParser::getBaseQuality()
{
    int *s = new int[source_bamread->core.l_qseq];
    int qual;
    for (int i = 0; i < source_bamread->core.l_qseq; ++i)
    {
        qual = bam_get_qual(source_bamread)[i];
        s[i] = qual;
        // std::cout << s[i] << std::endl;
    }

    return s;
}

uint8_t ReadParser::getMapQuality()
{
    return source_bamread->core.qual;
}

char complement(char n)
{
    switch (n)
    {
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    }

    return ' ';
}

std::string ReadParser::getReverseComplement(std::string nucs)
{
    transform(
        begin(nucs),
        end(nucs),
        begin(nucs),
        complement);

    std::reverse(nucs.begin(), nucs.end());

    return nucs;
}

void ReadParser::replaceToReverseComplement(std::string *nucs)
{
    transform(
        begin(*nucs),
        end(*nucs),
        begin(*nucs),
        complement);

    std::reverse(nucs->begin(), nucs->end());
}

int32_t ReadParser::getLengthSequence()
{
    return source_bamread->core.l_qseq;
}

bool ReadParser::isReverse()
{
    if ((source_bamread->core.flag & BAM_FREVERSE))
    {
        return true;
    }

    return false;
}

bool ReadParser::isMateReverse()
{
    if ((source_bamread->core.flag & BAM_FMREVERSE))
    {
        return true;
    }

    return false;
}

bool ReadParser::isFirstRead()
{
    if ((source_bamread->core.flag & BAM_FREAD1))
    {
        return true;
    }

    return false;
}

bool ReadParser::isSecondRead()
{
    if ((source_bamread->core.flag & BAM_FREAD2))
    {
        return true;
    }

    return false;
}

std::string ReadParser::getEdgeSeqFromStartSeq(uint32_t number)
{
    return getSequence().substr(0, number);
}

std::string ReadParser::getEdgeSeqFromEndSeq(uint32_t number)
{
    return getSequence().substr(source_bamread->core.l_qseq - number, source_bamread->core.l_qseq);
}

bool ReadParser::isNotPassingFilters()
{
    if ((source_bamread->core.flag & BAM_FQCFAIL))
    {
        return true;
    }

    return false;
}

bool ReadParser::isPCR()
{
    if ((source_bamread->core.flag & BAM_FDUP))
    {
        return true;
    }

    return false;
}
bool ReadParser::isSupplementaryAlignment()
{
    if ((source_bamread->core.flag & BAM_FSUPPLEMENTARY))
    {
        return true;
    }

    return false;
}

std::string ReadParser::getSoftClippedSequenceStart()
{
    auto cigar = getCigar();
    if (!(cigar.at(0).getOperatorName() == 'S'))
    {
        return "";
    }

    if (cigar.at(0).getLength() <= 2)
    {
        return "";
    }

    return getSequence().substr(0, cigar.at(0).getLength());
}

std::string ReadParser::getSoftClippedSequenceEnd()
{
    auto cigar = getCigar();
    if (!(cigar.at(cigar.size() - 1).getOperatorName() == 'S'))
    {
        return "";
    }

    if (cigar.at(cigar.size() - 1).getLength() <= 2)
    {
        return "";
    }

    return getSequence().substr(getSequence().size() - cigar.at(cigar.size() - 1).getLength(), getSequence().size());
}

int ReadParser::getLastToStartMissMatchPosMD()
{
    std::vector<ReadParser::AlignMD> md = getAlignMD();

    int sizeAcc = 0;
    // std::cout << "#md.size() :" << md.size() << std::endl;
    for (int i = md.size() - 1; i >= 0; i--)
    {
        AlignMD lastmd = md.at(i);
        sizeAcc += lastmd.size;
        // std::cout << "#" << lastmd.operate << " "<< lastmd.size << std::endl;
        if (lastmd.operate != 'M')
        {
            break;
        }
    }

    return sizeAcc;
}

int ReadParser::getSecondLastToStartMissMatchPosMD()
{
    std::vector<ReadParser::AlignMD> md = getAlignMD();

    int sizeAcc = 0;
    int numberFound = 0;
    for (int i = md.size() - 1; i >= 0; i--)
    {
        AlignMD lastmd = md.at(i);
        sizeAcc += lastmd.size;
        if (lastmd.operate != 'M')
        {
            numberFound++;
            if (numberFound >= 2)
            {
                break;
            }
        }
    }

    if (numberFound >= 2)
    {
        return sizeAcc;
    }

    return 0;
}

int ReadParser::getStartToEndMissMatchPosMD()
{
    std::vector<ReadParser::AlignMD> md = getAlignMD();

    int sizeAcc = 0;
    for (int i = 0; i < md.size(); i++)
    {
        AlignMD lastmd = md.at(i);
        sizeAcc += lastmd.size;
        // std::cout << "#" << lastmd.operate << " = " << lastmd.size << std::endl;
        if (lastmd.operate != 'M')
        {
            break;
        }
    }

    return sizeAcc;
}

int ReadParser::getSecondStartToEndMissMatchPosMD()
{
    std::vector<ReadParser::AlignMD> md = getAlignMD();

    int sizeAcc = 0;
    int numberFound = 0;
    for (int i = 0; i < md.size(); i++)
    {
        AlignMD lastmd = md.at(i);
        sizeAcc += lastmd.size;
        if (lastmd.operate != 'M')
        {
            numberFound++;
            if (numberFound >= 2)
            {
                break;
            }
        }
    }

    if (numberFound >= 2)
    {
        return sizeAcc;
    }

    return 0;
}

std::vector<ReadParser::AlignMD> ReadParser::getAlignMD()
{
    std::vector<ReadParser::AlignMD> alignMDs;
    const char *mdtagchar = "MD";
    auto aux = bam_aux_get(source_bamread, mdtagchar);
    auto auxChar = bam_aux2Z(aux);
    // std::cout << auxChar << std::endl;

    std::string collectNumber;
    std::string operate = "M";
    bool foundDel = false;
    int delCount = 0;
    for (int i = 0; i < strlen(auxChar); i++)
    {
        const char symbol = toupper(auxChar[i]);

        if (isdigit(symbol))
        {
            collectNumber += symbol;

            if (foundDel)
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'D';
                Aalignmd.size = delCount;
                alignMDs.push_back(Aalignmd);
                collectNumber = "";
                delCount = 0;
                foundDel = false;
            }
        }
        else
        {
            if (foundDel)
            {
                delCount++;
                continue;
            }

            if (collectNumber != "")
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'M';
                Aalignmd.size = std::stoi(collectNumber);
                alignMDs.push_back(Aalignmd);
                collectNumber = "";
            }

            if (symbol == '^')
            {
                foundDel = true;
                delCount = 0;
            }

            else if (symbol == 'A')
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'A';
                Aalignmd.size = 1;
                alignMDs.push_back(Aalignmd);
            }
            else if (symbol == 'G')
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'G';
                Aalignmd.size = 1;
                alignMDs.push_back(Aalignmd);
            }
            else if (symbol == 'T')
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'T';
                Aalignmd.size = 1;
                alignMDs.push_back(Aalignmd);
            }
            else if (symbol == 'C')
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'C';
                Aalignmd.size = 1;
                alignMDs.push_back(Aalignmd);
            }
        }
    }

    if (collectNumber != "")
    {
        ReadParser::AlignMD Aalignmd;
        Aalignmd.operate = 'M';
        Aalignmd.size = std::stoi(collectNumber);
        alignMDs.push_back(Aalignmd);
        collectNumber = "";
    }

    // for (auto n : alignMDs)
    // {
    //     std::cout << n.operate << " // " << n.size << std::endl;
    // }

    return alignMDs;
}

std::vector<ReadParser::SATag> ReadParser::getSATag()
{
    std::vector<SATag> saTag;

    // std::vector<ReadParser::AlignMD> alignMDs;
    const char *satagchar = "SA";
    auto aux = bam_aux_get(source_bamread, satagchar);
    if (aux == 0)
    {
        return saTag;
    }
    auto auxChar = bam_aux2Z(aux);
    std::string auxString(auxChar);

    std::vector<std::string> saString = splitText(auxString, ';');
    std::vector<std::string> subSplit;

    for (auto n : saString)
    {
        subSplit = splitText(n, ',');
        if (subSplit.size() != 0)
        {
            SATag tempSAtag;
            tempSAtag.chrname = subSplit.at(0);
            tempSAtag.pos = std::stol(subSplit.at(1), nullptr, 0);
            tempSAtag.strand = subSplit.at(2);
            tempSAtag.cigar = getCigarByString(subSplit.at(3));
            // std::cout << "cigar" << std::endl;
            // for (auto n:tempSAtag.cigar) {
            //     std::cout << n.getOperatorName() << " " << n.getLength() << std::endl;
            // }

            tempSAtag.mapQ = (uint8_t)atoi(subSplit.at(4).c_str());
            tempSAtag.NM = atoi(subSplit.at(5).c_str());

            saTag.push_back(tempSAtag);
        }
    }

    return saTag;
}

std::vector<std::string> ReadParser::splitText(std::string s, char delimiter)
{
    std::vector<std::string> tokens;

    std::stringstream ss;
    ss.str(s);

    std::string intermediate;

    while (getline(ss, intermediate, delimiter))
    {
        tokens.push_back(intermediate);
    }

    return tokens;
}

std::vector<ReadParser::Cigar> ReadParser::getCigarByString(std::string cigartext)
{
    std::vector<ReadParser::Cigar> cigar;
    std::string digitstring;
    for (auto n : cigartext)
    {
        if (isdigit(n))
        {
            digitstring += n;
        }
        else
        {
            Cigar tempCigar;
            tempCigar.setLength(atoi(digitstring.c_str()));
            digitstring = "";
            tempCigar.setOperatorName(n);
            cigar.push_back(tempCigar);
        }
    }

    return cigar;
}

bool ReadParser::isProperlyAligned()
{
    if ((source_bamread->core.flag & BAM_FPROPER_PAIR))
    {
        return true;
    }

    return false;
}