#include "evidence.h"
#include <cmath>

Evidence::Evidence()
{
    // mapqlist = new std::vector<uint8_t>;
}

Evidence::~Evidence()
{
    // associateReadLists.clear();
    // mapqlist->clear();
    // free (mapqlist);
}

std::string Evidence::getVariantType()
{
    return variantType;
}

void Evidence::setVariantType(std::string varianttype)
{
    variantType = varianttype;
}

void Evidence::setEndDiscordantRead(int32_t end)
{
    endDiscordantRead = end;
}

void Evidence::setPosDiscordantRead(int32_t pos)
{
    posDiscordantRead = pos;
}

int32_t Evidence::getPosDiscordantRead()
{
    return posDiscordantRead;
}

int32_t Evidence::getEndDiscordantRead()
{
    return endDiscordantRead;
}

std::string Evidence::getChr()
{
    return chr;
}

std::string Evidence::getEndChr()
{
    return mchr;
}

void Evidence::setFrequency(int count)
{
    frequency = count;
}

void Evidence::incrementFrequency()
{
    frequency++;
}

void Evidence::setChr(std::string t_chr)
{
    chr = t_chr;
}

void Evidence::setEndChr(std::string matechr)
{
    mchr = matechr;
}

void Evidence::setForwardDirection(bool isforward)
{
    forwardDirection = isforward;
}

bool Evidence::getForwardDirection()
{
    return forwardDirection;
}

void Evidence::setBackwardDirection(bool isbackward)
{
    backwardDirection = isbackward;
}

bool Evidence::getBackwardDirection()
{
    return backwardDirection;
}

int Evidence::getFrequency()
{
    return frequency;
}

void Evidence::addAssociateRead(int32_t pos, int32_t matepos)
{
    associateRead ar;
    ar.posDiscordantRead = pos;
    ar.mateposDiscordantRead = matepos;

    associateReadLists.push_back(ar);
}

//void Evidence::addErrorAssociateRead(int32_t pos, int32_t matepos)
//{
//    associateRead ar;
//    ar.posDiscordantRead = pos;
//    ar.mateposDiscordantRead = matepos;
//
//    errorAssociateReadLists.push_back(ar);
//}

void Evidence::setLastPosDiscordantRead(int32_t pos)
{
    lastposDiscordantRead = pos;
}

void Evidence::setLastEndDiscordantRead(int32_t end)
{
    lastendDiscordantRead = end;
}

int32_t Evidence::getLastPosDiscordantRead()
{
    return lastposDiscordantRead;
}

int32_t Evidence::getLastEndDiscordantRead()
{
    return lastendDiscordantRead;
}

void Evidence::addMapQ(uint8_t mapq)
{
    (mapqlist).push_back(mapq);
}

std::vector<uint8_t> *Evidence::getMapQVector()
{
    return &mapqlist;
}

std::string Evidence::convertMapQlistToCommaString(std::vector<uint8_t> *mapqlist)
{
    std::string mapqString;

    bool first = false;
    for (auto n : *mapqlist)
    {
        if (!first)
        {
            mapqString.append(std::to_string(n));
            first = true;
        }
        else
        {
            mapqString.append("," + std::to_string(n));
        }
    }

    return mapqString;
}

uint8_t Evidence::getMaxMapQ()
{
    uint8_t max = 0;
    for (auto n : mapqlist)
    {
        if (n > max)
        {
            max = n;
        }
    }

    return max;
}

uint8_t Evidence::getMinMapQ()
{
    uint8_t min = 255;
    for (auto n : mapqlist)
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

uint8_t Evidence::getAvgMapQ()
{
    uint8_t countQual = 0;
    for (auto n : mapqlist)
    {
        countQual += n;
    }

    return countQual / mapqlist.size();
}

uint8_t Evidence::getMaxRPMapQ()
{
    uint8_t max = 0;
    for (auto n : rpmapqlist)
    {
        if (n > max)
        {
            max = n;
        }
    }

    return max;
}

uint8_t Evidence::getMinRPMapQ()
{
    uint8_t min = 255;
    for (auto n : rpmapqlist)
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

int Evidence::getNumberOfRP()
{
    return rpmapqlist.size();
}

uint8_t Evidence::getAvgRPMapQ()
{
    uint8_t countQual = 0;
    for (auto n : rpmapqlist)
    {
        countQual += n;
    }

    return countQual / rpmapqlist.size();
}

std::vector<std::string> split(const std::string &s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

void Evidence::setEvidenceByString(std::string line)
{
    // std::cout << "setEvidenceByString" << std::endl;
    std::vector<std::string> results = split(line, '\t');
    int i = 0;
    Evidence e;
    int32_t pv_pos = 0;
    for (auto n : results)
    {
        if (i == 0)
        {
            //CHROM
            setChr(n);
            setEndChr(n);
        }
        else if (i == 1)
        {
            //POS
            pv_pos = std::stol(n, nullptr, 0);
            setPos(pv_pos);
            // setPosDiscordantRead(ul);
            // setChr(n);
            // std::cout << n << "/" << getChr() << std::endl;
        }
        else if (i == 2)
        {
            //ID
            // setID
            // uint32_t ul = std::stol(n, nullptr, 0);
            // setPosDiscordantRead(ul);
            // std::cout << n << "/" << getPosDiscordantRead() << std::endl;
        }
        else if (i == 3)
        {
            //REF
            // uint32_t ul = std::stol(n, nullptr, 0);
            // setLastPosDiscordantRead(ul);
            // std::cout << n << "/" << getLastPosDiscordantRead() << std::endl;
        }
        else if (i == 4)
        {
            //ALT
            // setEndChr(n);
        }
        else if (i == 5)
        {
            //QUAL
            // uint32_t ul = std::stol(n, nullptr, 0);
            // setEndDiscordantRead(ul);
            // std::cout << n << "/" << getEndDiscordantRead() << std::endl;
        }
        else if (i == 6)
        {
            //FILTER
            // uint32_t ul = std::stol(n, nullptr, 0);
            // setLastEndDiscordantRead(ul);
            // std::cout << n << "/" << getLastEndDiscordantRead() << std::endl;
        }
        else if (i == 7)
        {
            //INFO
            std::vector<std::string> infolist = split(n, ';');
            //            int32_t negativeCIPOS = 0;
            //            int32_t postitiveCIPOS = 0;
            //            int32_t negativeCIEND = 0;
            //            int32_t postitiveCIEND = 0;
            int32_t pv_end = 0;

            for (std::string ainfo : infolist)
            {
                // std::cout << ainfo << "  -  " << getKeybyText(ainfo)  << std::endl;

                if (getKeybyText(ainfo) == "END")
                {
                    pv_end = std::stol(getValuebyText(ainfo), nullptr, 0);
                    setEnd(pv_end);
                    continue;
                }
                else if (getKeybyText(ainfo) == "CIPOS")
                {

                    std::vector<std::string> nineSection = split(getValuebyText(ainfo), ',');
                    if (nineSection.size() == 2)
                    {
                        //                        negativeCIPOS = std::stol(nineSection[0], nullptr, 0);
                        //                        postitiveCIPOS = std::stol(nineSection[1], nullptr, 0);
                        setCiPosLeft(std::stol(nineSection[0], nullptr, 0));
                        setCiPosRight(std::stol(nineSection[1], nullptr, 0));
                    }
                    continue;
                }
                else if (getKeybyText(ainfo) == "CIEND")
                {
                    std::string pv_CIEND = getValuebyText(ainfo);
                    std::vector<std::string> nineSection = split(pv_CIEND, ',');
                    if (nineSection.size() == 2)
                    {
                        //                        negativeCIEND = std::stol(nineSection[0], nullptr, 0);
                        //                        postitiveCIEND = std::stol(nineSection[1], nullptr, 0);
                        setCiEndLeft(std::stol(nineSection[0], nullptr, 0));
                        setCiEndRight(std::stol(nineSection[1], nullptr, 0));
                    }
                    continue;
                }
                else if (getKeybyText(ainfo) == "SVTYPE")
                {
                    setVariantType(getValuebyText(ainfo));
                    continue;
                }
                else if (getKeybyText(ainfo) == "CHR2")
                {
                    setEndChr(getValuebyText(ainfo));
                    continue;
                }
                else if (getKeybyText(ainfo) == "iPRIns_MQL")
                {
                    std::string pv_BMQL = getValuebyText(ainfo);
                    std::vector<std::string> nineSection = split(getValuebyText(ainfo), ',');
                    for (std::string a_mq : nineSection)
                    {
                        uint8_t t_mq = (uint8_t)atoi(a_mq.c_str());
                        addMapQ(t_mq);
                    }
                    continue;
                }
                else if (getKeybyText(ainfo) == "iPRIns_RPMQL")
                {
                    std::string pv_BMQL = getValuebyText(ainfo);
                    std::vector<std::string> nineSection = split(getValuebyText(ainfo), ',');
                    for (std::string a_mq : nineSection)
                    {
                        uint8_t t_mq = (uint8_t)atoi(a_mq.c_str());
                        rpmapqlist.push_back(t_mq);
                    }
                    continue;
                }
                else if (getKeybyText(ainfo) == "iPRIns_POS_DR")
                {
                    std::string pv_BMQL = getValuebyText(ainfo);
                    std::vector<std::string> nineSection = split(getValuebyText(ainfo), ',');
                    setPosDiscordantRead(std::stol(nineSection.at(0), nullptr, 0));
                    setLastPosDiscordantRead(std::stol(nineSection.at(1), nullptr, 0));
                    continue;
                }
                else if (getKeybyText(ainfo) == "iPRIns_END_DR")
                {
                    std::string pv_BMQL = getValuebyText(ainfo);
                    std::vector<std::string> nineSection = split(getValuebyText(ainfo), ',');
                    setEndDiscordantRead(std::stol(nineSection.at(0), nullptr, 0));
                    setLastEndDiscordantRead(std::stol(nineSection.at(1), nullptr, 0));
                    continue;
                }
 
                else if (getKeybyText(ainfo) == "iPRIns_LNGMATCH")
                {
                    LNGMATCH = atoi(getValuebyText(ainfo).c_str());
                    continue;
                }

                else if (getKeybyText(ainfo) == "iPRIns_FREQ")
                {
                    setFrequency(atoi(getValuebyText(ainfo).c_str()));
                    continue;
                }

                else if (getKeybyText(ainfo) == "iPRIns_MARK")
                {
                    setMark(getValuebyText(ainfo));
                    continue;
                }

                else if (getKeybyText(ainfo) == "iPRIns_SA")
                {
                    std::vector<std::string> nineSection = split(getValuebyText(ainfo), ',');
                    for (auto m : nineSection)
                    {
                        std::vector<std::string> n = split(m, ':');
                        // for (auto n : BAltSA)
                        // {
                        addAlterSA(n.at(0), std::stol(n.at(1), nullptr, 0));
                        // }
                    }
                    continue;
                }

                else if (getKeybyText(ainfo) == "iPRIns_EPOS")
                {
                    evidencefrom = std::stol(getValuebyText(ainfo), nullptr, 0);
                    continue;
                }

            }
        }

        i++;
    }
}

std::string Evidence::getKeybyText(std::string n)
{
    std::vector<std::string> nlist = split(n, '=');

    if (nlist.size() == 2)
    {
        return nlist[0];
    }

    return "";
}

std::string Evidence::getID()
{
    return ID;
}

std::string Evidence::getRef()
{
    return ref;
}

std::string Evidence::getAlt()
{
    return alt;
}

uint8_t Evidence::getQual()
{
    return qual;
}

std::string Evidence::getFilter()
{
    return filter;
}

std::string Evidence::getInfoString()
{
    std::string result;
    result.append("END=" + std::to_string(getEnd()) + ";");
    result.append("SVTYPE=" + getVariantType() + ";");

    if (getVariantType() != "BND")
    {
        result.append("SVLEN=" + std::to_string(getEnd() - getPos()) + ";");
    }
    if (getVariantType() == "BND")
    {
        result.append("CHR2=" + getEndChr() + ";");
    }

    // if (getComment() != "")
    // {
    //     result.append("COMMENT=" + getComment() + ";");
    // }
    if (getFrequency() != 0)
    {
        result.append("iPRIns_FREQ=" + std::to_string(getFrequency()) + ";");
    }
    result.append("iPRIns_LNGMATCH=" + std::to_string(LNGMATCH) + ";");
    result.append("iPRIns_MQL=" + convertMapQlistToCommaString(getMapQVector()) + ";");
    if (getMark() != "")
    {
        result.append("iPRIns_MARK=" + getMark() + ";");
    }
    if (getRPMapQ()->size() != 0)
    {
        result.append("iPRIns_RPMQL=" + convertMapQlistToCommaString(getRPMapQ()) + ";");
    }

    if (getCiPosLeft() != 0 || getCiPosRight() != 0)
    {
        result.append("CIPOS=" + std::to_string(getCiPosLeft()) + "," + std::to_string(getCiPosRight()) + ";");
    }
    if (getCiEndLeft() != 0 || getCiEndRight() != 0)
    {
        result.append("CIEND=" + std::to_string(getCiEndLeft()) + "," + std::to_string(getCiEndRight()) + ";");
    }

    if (getPosDiscordantRead() != 0 || getLastPosDiscordantRead() != 0)
    {
        result.append("iPRIns_POS_DR=" + std::to_string(getPosDiscordantRead()) + "," + std::to_string(getLastPosDiscordantRead()) + ";");
    }

    if (getEndDiscordantRead() != 0 || getLastEndDiscordantRead() != 0)
    {
        result.append("iPRIns_END_DR=" + std::to_string(getEndDiscordantRead()) + "," + std::to_string(getLastEndDiscordantRead()) + ";");
    }

    if (AltSA.size() != 0)
    {
        result.append("iPRIns_SA=");
        bool first = true;
        for (auto n : AltSA)
        {
            if (first)
            {
                first = false;
                result.append(n.chr + ":" + std::to_string(n.pos));
            }
            result.append("," + n.chr + ":" + std::to_string(n.pos));
        }
        result.append(";");
    }

    if (evidencefrom != 0)
    {
        result.append("iPRIns_EPOS=" + std::to_string(evidencefrom) + ";");
    }

    return result;
}

void Evidence::addAlterSA(std::string chr, int32_t pos)
{
    AlternativeSA tempALT;
    tempALT.chr = chr;
    tempALT.pos = pos;
    AltSA.push_back(tempALT);
}
void Evidence::setEvidenceFrom(int32_t pos)
{
    evidencefrom = pos;
}

int32_t Evidence::getEvidencePos()
{
    return evidencefrom;
}

int Evidence::countEvidencePosNear(int32_t pos, int32_t overlapped)
{
    int count = 0;
    for (auto n : AltSA)
    {
        if (checkBetween(n.pos, pos, overlapped))
        {
            count++;
        }
    }

    return count;
}

int Evidence::countDiffEvidencePos(int32_t overlapped)
{
    std::vector<AlternativeSA> tempAltSA;
    for (auto n : AltSA)
    {
        bool found = false;
        for (auto m : tempAltSA)
        {
            if (m.chr != n.chr)
            {
                continue;
            }

            if (checkBetween(n.pos, m.pos, overlapped))
            {
                found = true;
                break;
            }
        }

        if (!found) 
        {
            tempAltSA.push_back(n);
        }
    }

    return tempAltSA.size();
}

int Evidence::getAltSASize()
{
    return AltSA.size();
}

bool Evidence::checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped)
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

void Evidence::setMark(std::string mark)
{
    Evidence::mark = mark;
}

std::string Evidence::getMark()
{
    return mark;
}

// std::string Evidence::getComment()
// {
//     return comment;
// }

// void Evidence::setComment(std::string comment)
// {
//     Evidence::comment = comment;
// }

void Evidence::setMapQList(std::vector<uint8_t> mapqs)
{
    mapqlist = mapqs;
}

void Evidence::setRPMapQ(std::vector<uint8_t> rpmapq)
{
    Evidence::rpmapqlist = rpmapq;
}

std::vector<uint8_t> *Evidence::getRPMapQ()
{
    return &rpmapqlist;
}

void Evidence::setID(std::string id)
{
    ID = id;
}

std::string Evidence::getResultVcfFormatString()
{
    std::string result;

    // 1. #CHROM
    result.append(getChr());

    // 2. POS
    result.append("\t");
    result.append(std::to_string(getPos()));

    // 3. ID
    result.append("\t");
    if (getID() == "")
    {
        result.append(".");
    }
    else
    {
        result.append(getID());
    }

    // 4. REF
    result.append("\t");
    if (getRef() == "")
    {
        result.append(".");
    }
    else
    {
        result.append(getRef());
    }

    // 5. ALT
    result.append("\t");
    if (getAlt() == "")
    {
        result.append(".");
    }
    else
    {
        result.append(getAlt());
    }

    // 6. QUAL
    result.append("\t");
    result.append(std::to_string(getQual()));

    // 7. FILTER
    result.append("\t");
    if (getFilter() == "")
    {
        result.append(".");
    }
    else
    {
        result.append(getFilter());
    }

    // 8. INFO
    result.append("\t");
    result.append(getInfoString());

    return result;
}

std::string Evidence::getValuebyText(std::string n)
{
    std::vector<std::string> nlist = split(n, '=');

    if (nlist.size() == 2)
    {
        return nlist[1];
    }

    return "";
}

std::string Evidence::convertToVcfString()
{
    std::string buf;
    //CHROM
    buf.append(getChr() + "\t");
    //POS
    buf.append(std::to_string(getPos()) + "\t");
    //ID
    buf.append(".\t");
    //REF
    buf.append(".\t");
    //ALT
    buf.append(".\t");
    //QUAL
    buf.append(".\t");
    //FILTER
    buf.append(".\t");
    //INFO
    buf.append("END=" + std::to_string(getEnd()) + ";");
    buf.append("SVTYPE=" + getVariantType() + ";");
    buf.append("SVLEN=" + std::to_string(getSvLength()) + ";");
    buf.append("CIPOS=" + std::to_string(getCiPosLeft()) + "," + std::to_string(getCiPosRight()) + ";");
    buf.append("CIEND=" + std::to_string(getCiEndLeft()) + "," + std::to_string(getCiEndRight()) + ";");
    buf.append("CHR2=" + getEndChr() + ";");
    buf.append("iPRIns_MQL=" + convertMapQlistToCommaString(getMapQVector()) + ";");
    buf.append("iPRIns_POS_DR=" + std::to_string(getPosDiscordantRead()) + "," + std::to_string(getLastPosDiscordantRead()) + ";");
    buf.append("iPRIns_END_DR=" + std::to_string(getEndDiscordantRead()) + "," + std::to_string(getLastEndDiscordantRead()) + ";");
    if (getMark() != "")
    {
        buf.append("iPRIns_MARK=" + getMark() + ";");
    }
    if (getFrequency() != 0)
    {
        buf.append("iPRIns_FREQ=" + std::to_string(getFrequency()) + ";");
    }
    if (rpmapqlist.size() != 0)
    {
        buf.append("iPRIns_RPMQL=" + convertMapQlistToCommaString(&rpmapqlist) + ";");
    }

    return buf;
}

int32_t Evidence::getSvLength()
{
    return getEnd() - getPos();
}

void Evidence::setQuailtyPass(bool QuailtyPass)
{
    Evidence::QuailtyPass = QuailtyPass;
}

bool Evidence::isQuailtyPass() const
{
    return QuailtyPass;
}

std::string Evidence::getSVType()
{
    return variantType;
}

int32_t Evidence::getPos()
{
    return pos;
}

void Evidence::setPos(int32_t pos)
{
    Evidence::pos = pos;
}

int32_t Evidence::getEnd()
{
    return end;
}

void Evidence::setEnd(int32_t end)
{
    Evidence::end = end;
}

int32_t Evidence::getLastPos()
{
    return lastpos;
}

void Evidence::setLastPos(int32_t lastpos)
{
    Evidence::lastpos = lastpos;
}

int32_t Evidence::getLastEnd()
{
    return lastend;
}

void Evidence::setLastEnd(int32_t lastend)
{
    Evidence::lastend = lastend;
}

int32_t Evidence::getCiPosLeft()
{

    return ciPosLeft;
}

void Evidence::setCiPosLeft(int32_t ciPosLeft)
{
    Evidence::ciPosLeft = ciPosLeft;
}

int32_t Evidence::getCiPosRight()
{
    return ciPosRight;
}

void Evidence::setCiPosRight(int32_t ciPosRight)
{
    Evidence::ciPosRight = ciPosRight;
}

int32_t Evidence::getCiEndLeft()
{
    return ciEndLeft;
}

void Evidence::setCiEndLeft(int32_t ciEndLeft)
{
    Evidence::ciEndLeft = ciEndLeft;
}

int32_t Evidence::getCiEndRight()
{
    return ciEndRight;
}

void Evidence::setCiEndRight(int32_t ciEndRight)
{
    Evidence::ciEndRight = ciEndRight;
}

std::vector<std::string> &Evidence::getMultipleEndChromosome()
{
    return multipleEndChromosome;
}

void Evidence::setMultipleEndChromosome(const std::vector<std::string> &multipleEndChromosome)
{
    Evidence::multipleEndChromosome = multipleEndChromosome;
}

int Evidence::getNumberOfZeroMapQ()
{
    int count = 0;
    for (auto n : mapqlist)
    {
        if (n == 0)
        {
            count++;
        }
    }

    return count;
}


int Evidence::getNumberOfZeroRPMapQ()
{
    int count = 0;
    for (auto n : rpmapqlist)
    {
        if (n == 0)
        {
            count++;
        }
    }

    return count;
}

bool Evidence::haveSomeMapQMoreThan(uint8_t qual)
{
    for (auto n : mapqlist)
    {
        if (n > qual)
        {
            return true;
        }
    }

    return false;
}

bool Evidence::haveSomeMapQLessThan(uint8_t qual)
{
    for (auto n : mapqlist)
    {
        if (n < qual)
        {
            return true;
        }
    }

    return false;
}
