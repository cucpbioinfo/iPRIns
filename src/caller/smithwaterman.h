#ifndef SMITHWATERMAN_H
#define SMITHWATERMAN_H

#include <cstdint>
#include <string>
#include <vector>

class SmithWaterman
{

  public:

    struct ScoreAlignment
    {
      int32_t pos=0;
      int32_t end=0;
      int32_t posseq = 0;
      int32_t endseq = 0;
      std::string pattern; // Pattern of sequence read
      int lastscore=0;
      int scorepattern = 0;
      int countpatterndel = 0;
      int countpatternins = 0;
      int countpatternmatch = 0;
      int countpatternmissmatch = 0;
    };

    SmithWaterman();
    SmithWaterman(std::string *ref, int32_t position,bool startToEnd);

    void RelocateSequenceDeletionFirst(std::string sequence, int32_t startposition);
    void RelocateSequenceDeletionLast();
    int calculateScorePattern(ScoreAlignment *scorealignment);
    void countScorePattern(ScoreAlignment *scorealignment);

    std::vector<ScoreAlignment> findEndToStart(std::string seq);
    std::vector<ScoreAlignment> findStartToEnd(std::string seq);

    bool refineScoreAlignment(SmithWaterman::ScoreAlignment *scorealignment,int numberofalign);
    void reverseReference();
    void getResult();

    int32_t getRefPos();

    void setRefPos(int32_t RefPos);

    std::string getReference();

    void setReference(std::string *reference);
    std::vector<ScoreAlignment> removeRedundant(std::vector<ScoreAlignment> *scorelist);

    int findMaxMatchInsertion(std::string *query);

  private:
    std::string reference;
    int32_t RefPos;

    void makeScoreMatrix(std::vector<std::vector<int>> *cell,std::string *seq);
    void makeScoreRow(std::vector<int> *previousCell,std::vector<int> *currentCell,char letter);

    int matchScore = 1;
    int missmatchScore = -1;
    int gapPenalty  = -5;
    int getMaxInteger(std::vector<int> *score);
    ScoreAlignment getNextPath(std::vector<std::vector<int>> *cell,std::string *seq,int i,int j,bool findStartToEnd);

    int getScoreUpper(std::vector<std::vector<int>> *cell,int i,int j);
    int getScoreUpperleft(std::vector<std::vector<int>> *cell,int i,int j);
    int getScoreLeft(std::vector<std::vector<int>> *cell,int i,int j);
    void goToUpper(int *i,int *j);
    void goToUpperleft(int *i,int *j);
    void goToLeft(int *i,int *j);

    SmithWaterman::ScoreAlignment calculatePositionfindEndToStart(int pos,std::string pattern);
    SmithWaterman::ScoreAlignment calculatePositionfindStartToEnd(int pos,std::string pattern);

};

#endif
