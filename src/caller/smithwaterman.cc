#include "smithwaterman.h"
#include <iostream>
#include <algorithm>


SmithWaterman::SmithWaterman()
{
}

SmithWaterman::SmithWaterman(std::string *ref, int32_t position, bool startToEnd)
{
    setReference(ref);
    setRefPos(position);
    if (startToEnd)
    {
        reverseReference();
    }
}

std::string SmithWaterman::getReference()
{
    return reference;
}

void SmithWaterman::reverseReference()
{
    reverse(reference.begin(), reference.end());
}

void SmithWaterman::setReference(std::string *reference)
{
    SmithWaterman::reference = *reference;
}

int32_t SmithWaterman::getRefPos()
{
    return RefPos;
}

void SmithWaterman::setRefPos(int32_t RefPos)
{
    SmithWaterman::RefPos = RefPos;
}

int SmithWaterman::findMaxMatchInsertion(std::string *seq)
{
    std::vector<std::vector<int>> cell;
    makeScoreMatrix(&cell, seq);

    int maxMatch=0;
    for (std::vector<int> n : cell)
    {
        for (int m : n)
        {
            if (m>maxMatch)
            {
                maxMatch = m;
            }
            
        }
    }


    return maxMatch;
}

std::vector<SmithWaterman::ScoreAlignment> SmithWaterman::findEndToStart(std::string seq)
{
    std::vector<std::vector<int>> cell;
    makeScoreMatrix(&cell, &seq);
    std::vector<SmithWaterman::ScoreAlignment> result;

    for (int i = 0; i < reference.length() + 1; i++)
    {
        int score = cell.at(seq.length()).at(i);
        if (score >= 4)
        {
            ScoreAlignment aResult = getNextPath(&cell, &seq, seq.size(), i, false);
            // countScorePattern(&aResult);
            // if (aResult.countpatternmatch >= 2)
            // {
            // aResult.pos = aResult.pos;
            // aResult.end = aResult.end;
            // aResult.scorepattern = calculateScorePattern(&aResult);
            // aResult.pa
            // std::cout << getRefPos() << std::endl;
            // std::cout << aResult.pattern << std::endl;

            aResult.lastscore = score;
            aResult.posseq = seq.size() - aResult.pattern.size();
            aResult.endseq = seq.size();

            result.push_back(aResult);
            // }
        }
    }

    return removeRedundant(&result);
}

std::vector<SmithWaterman::ScoreAlignment> SmithWaterman::findStartToEnd(std::string seq)
{
    reverse(seq.begin(), seq.end());
    std::vector<std::vector<int>> cell;
    makeScoreMatrix(&cell, &seq);
    std::vector<SmithWaterman::ScoreAlignment> result;

    for (int i = 0; i < reference.length() + 1; i++)
    {
        int score = cell.at(seq.length()).at(i);
        if (score >= 4)
        {
            ScoreAlignment aResult = getNextPath(&cell, &seq, seq.size(), i, true);
            // countScorePattern(&aResult);
            // if (aResult.countpatternmatch >= 3)
            // {

            // aResult.scorepattern = calculateScorePattern(&aResult);
            // int32_t end = getReference().size() - aResult.pos;
            // int32_t pos = getReference().size() - aResult.end;
            // aResult.pos = pos;
            // aResult.end = end;
            //  std::cout << "aResult.pos : " << aResult.pos << std::endl;
            // std::cout << "aResult.end  : " << aResult.end  << std::endl;

            // std::cout << "aResult.scorepattern : " << aResult.scorepattern << std::endl;
            // std::cout << "aResult.pattern : " << aResult.pattern << std::endl;
            // aResult.lastscore = score;

            // aResult.posseq = 0;
            // aResult.endseq = aResult.posseq + aResult.pattern.size();
            result.push_back(aResult);
            // }
        }
    }

    return removeRedundant(&result);
}

bool SmithWaterman::refineScoreAlignment(SmithWaterman::ScoreAlignment *scorealignment, int numberofalign)
{
    // return true;
    // if (numberofalign < scorealignment->pattern.size())
    // {
    //     int diff = scorealignment->pattern.size() - numberofalign;
    //     scorealignment->pattern = scorealignment->pattern.substr(0, numberofalign);
    //     scorealignment->endseq = scorealignment->endseq - diff;
    //     scorealignment->end = scorealignment->end - diff;
    //     return true;
    // }
    // else if (numberofalign == scorealignment->pattern.size())
    // {
    //     return true;
    // }

    if (numberofalign == scorealignment->pattern.size())
    {
        return true;
    }

    return false;
}

std::vector<SmithWaterman::ScoreAlignment> SmithWaterman::removeRedundant(std::vector<SmithWaterman::ScoreAlignment> *scorelist)
{
    // std::vector<ScoreAlignment> result;
    // for (auto n : *scorelist)
    // {
    //     bool added = false;
    //     for (auto m : result)
    //     {
    //         if (m.pos == n.pos && m.end == n.end)
    //         {
    //             added = true;
    //             break;
    //         }
    //     }

    //     if (!added)
    //     {
    //         result.push_back(n);
    //     }
    // }

    return *scorelist;
}

int SmithWaterman::calculateScorePattern(ScoreAlignment *ScoreAlignment)
{
    int score = 0;
    for (auto n : ScoreAlignment->pattern)
    {
        switch (n)
        {
        case 'M':
            score += 1;
            break;
        case 'X':
            score -= 1;
            break;
        case 'D':
            score -= 2;
            break;
        case 'I':
            score -= 2;
            break;
        default:
            break;
        }
    }

    return score;
}

void SmithWaterman::countScorePattern(ScoreAlignment *ScoreAlignment)
{

    for (auto n : ScoreAlignment->pattern)
    {
        switch (n)
        {
        case 'M':
            ScoreAlignment->countpatternmatch += 1;
            break;
        case 'X':
            ScoreAlignment->countpatternmissmatch += 1;
            break;
        case 'D':
            ScoreAlignment->countpatterndel += 1;
            break;
        case 'I':
            ScoreAlignment->countpatternins += 1;
            break;
        default:
            break;
        }
    }
}

int SmithWaterman::getScoreUpper(std::vector<std::vector<int>> *cell, int i, int j)
{
    return cell->at(i).at(j - 1);
}
int SmithWaterman::getScoreUpperleft(std::vector<std::vector<int>> *cell, int i, int j)
{
    return cell->at(i - 1).at(j - 1);
}
int SmithWaterman::getScoreLeft(std::vector<std::vector<int>> *cell, int i, int j)
{
    return cell->at(i - 1).at(j);
}

void SmithWaterman::goToUpper(int *i, int *j)
{
    *j = *j - 1;
}
void SmithWaterman::goToUpperleft(int *i, int *j)
{
    *i = *i - 1;
    *j = *j - 1;
}
void SmithWaterman::goToLeft(int *i, int *j)
{
    *i = *i - 1;
}

SmithWaterman::ScoreAlignment SmithWaterman::getNextPath(std::vector<std::vector<int>> *cell, std::string *seq, int i, int j, bool findStartToEnd)
{
    std::string pattern;
    int32_t start = j;

    int currentpos = cell->at(i).at(j);
    if (seq->at(i - 1) == reference.at(j - 1))
    {
        pattern += "M";
    }
    else
    {
        pattern += "X";
    }
    // std::cout << "seq : " << seq->at(i - 1) << " ref : " << reference.at(j - 1) << std::endl;
    // goToUpperleft(&i, &j);

    // if (currentpos != 9)
    // {
    //     ScoreAlignment score;
    //     return score;
    // }

    // std::cout << "i : " << i << " j : " << j << std::endl;
    int stateMissmatch = 0;
    int missmatch = 0;

    for (;;)
    {
        // std::cout << "----start --- " << std::endl;

        if (i == 1 || j == 1)
        {
            break;
        }

        if (stateMissmatch > 3)
        {
            break;
        }

        if (missmatch == 2)
        {
            break;
        }

        int currentpos = cell->at(i).at(j);

        // std::cout << "currentpos : " << currentpos << " current : " << seq->at(i-1) << std::endl;
        // std::cout << "current seq : " << seq->at(i - 1) << " ref : " << reference.at(j - 1) << std::endl;
        // std::cout << "current pattern : " << pattern << std::endl;

        int upper = getScoreUpper(cell, i, j);
        int upperleft = getScoreUpperleft(cell, i, j);
        int left = getScoreLeft(cell, i, j);

        // std::cout << " > upper : " << upper << std::endl;
        // std::cout << " > upperleft : " << upperleft << std::endl;
        // std::cout << " > left : " << left << std::endl;

        if (upper == 0 && upperleft == 0 && left == 0)
        {
            break;
        }

        if (upperleft != 0)
        {
            // if (upperleft == currentpos - 1)
            // {
            // std::cout << "upperleft" << std::endl;
            if (seq->at(i - 1) == reference.at(j - 1))
            {
                pattern += "M";
                stateMissmatch = 0;

                goToUpperleft(&i, &j);
                // std::cout << " > goto : goToUpperleft M" << std::endl;
                continue;
            }
            // }

            if (upperleft > upper && upperleft > left)
            {

                pattern += "X";
                stateMissmatch++;
                missmatch++;
                goToUpperleft(&i, &j);
                // std::cout << " > goto : goToUpperleft X" << std::endl;
                continue;
            }
        }

        if (upper > left)
        {
            pattern += "D";
            stateMissmatch++;
            missmatch++;
            goToUpper(&i, &j);
            //  std::cout << " > goto : goToUpper" << std::endl;
            continue;
        }

        pattern += "I";
        stateMissmatch++;
        missmatch++;
        goToLeft(&i, &j);
        // std::cout << " > goto : goToLeft" << std::endl;
        continue;
    }

    if (findStartToEnd)
    {
        return calculatePositionfindEndToStart(start, pattern);
    }
    else
    {
        return calculatePositionfindStartToEnd(start, pattern);
    }
}

SmithWaterman::ScoreAlignment SmithWaterman::calculatePositionfindStartToEnd(int pos, std::string pattern)
{

    ScoreAlignment score;

    int refRun = 0;
    for (auto n : pattern)
    {
        if (n == 'M')
        {
            refRun++;
        }
        else if (n == 'X')
        {
            refRun++;
        }
        else if (n == 'D')
        {
            refRun++;
        }
        else if (n == 'I')
        {
        }
    }

    int seqRun = 0;
    for (auto n : pattern)
    {
        if (n == 'M')
        {
            seqRun++;
        }
        else if (n == 'X')
        {
            seqRun++;
        }
        else if (n == 'D')
        {
        }
        else if (n == 'I')
        {
            seqRun++;
        }
    }

    score.endseq = seqRun;
    score.pos = pos - refRun;
    score.end = score.pos + refRun;
    score.pattern = pattern;

    // std::cout << "### pos : " << pos << " pattern : " << pattern << " refRun : " << refRun << " ref size : " << getReference().size() << std::endl;

    return score;
}

SmithWaterman::ScoreAlignment SmithWaterman::calculatePositionfindEndToStart(int pos, std::string pattern)
{

    ScoreAlignment score;

    int refRun = 0;
    for (auto n : pattern)
    {
        if (n == 'M')
        {
            refRun++;
        }
        else if (n == 'X')
        {
            refRun++;
        }
        else if (n == 'D')
        {
            refRun++;
        }
        else if (n == 'I')
        {
        }
    }

    int seqRun = 0;
    for (auto n : pattern)
    {
        if (n == 'M')
        {
            seqRun++;
        }
        else if (n == 'X')
        {
            seqRun++;
        }
        else if (n == 'D')
        {
        }
        else if (n == 'I')
        {
            seqRun++;
        }
    }

    score.endseq = seqRun;
    score.pos = getReference().size() - pos;
    score.end = score.pos + refRun;
    score.pattern = pattern;

    // std::cout << "### pos : " << pos << " pattern : " << pattern << " refRun : " << refRun << " ref size : " << getReference().size() << std::endl;

    return score;
}

void SmithWaterman::makeScoreMatrix(std::vector<std::vector<int>> *cell, std::string *seq)
{

    // std::cout << *seq << std::endl;
    // std::cout << reference.size() << std::endl;
    // make first row
    std::vector<int> firstRow;
    firstRow.push_back(0);
    // std::cout << "  # ";
    for (int i = 0; i < reference.size(); i++)
    {
        // std::cout << reference.at(i) << " ";
        firstRow.push_back(0);
    }
    cell->push_back(firstRow);

    // inner data
    std::vector<int> secondRow;

    // std::cout << std::endl;
    // std::cout << "  ";
    // for (auto n : firstRow)
    // {
    //     std::cout << n << " ";
    // }
    // std::cout << std::endl;

    for (auto letter : *seq)
    {
        makeScoreRow(&cell->back(), &secondRow, letter);
        cell->push_back(secondRow);
        // std::cout << letter << " ";
        // for (auto n : secondRow)
        // {
        //     std::cout << n << " ";
        // }

        // std::cout << std::endl;

        secondRow.clear();
    }

    // std::cout << std::endl;
}

void SmithWaterman::makeScoreRow(std::vector<int> *previousCell, std::vector<int> *currentCell, char letter)
{
    // std::cout << "---" << std::endl;
    int currentPreviousScore = 0;
    currentCell->push_back(currentPreviousScore);
    // std::cout << currentPreviousScore << " ";

    std::vector<int> tempScore;
    int selfScore = missmatchScore;
    for (int i = 0; i < reference.size(); i++)
    {
        // self score
        if (reference.at(i) == letter)
        {
            selfScore = matchScore;
        }
        tempScore.push_back(selfScore);
        // std::cout << reference.at(i) << " = " << letter << std::endl;
        // std::cout << std::endl << reference.at(i) << " = " << letter << " " << selfScore << "  :  ";
        // previous score
        tempScore.push_back(currentPreviousScore + gapPenalty);

        // upper left
        tempScore.push_back(previousCell->at(i) + selfScore);

        // upper
        tempScore.push_back(previousCell->at(i + 1) + gapPenalty);

        currentPreviousScore = getMaxInteger(&tempScore);
        currentCell->push_back(currentPreviousScore);

        // std::cout << currentPreviousScore << " ";
        //reset variable
        selfScore = missmatchScore;
        tempScore.clear();
    }
}

int SmithWaterman::getMaxInteger(std::vector<int> *score)
{
    int max = 0;
    for (auto n : *score)
    {
        if (n > max)
        {
            max = n;
        }
    }

    return max;
}
