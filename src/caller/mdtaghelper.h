#ifndef MDTAGHELPER_H
#define MDTAGHELPER_H

#include <string>
#include <vector>

class MDtagHelper {
private:

public:
    MDtagHelper(std::string mdtext,std::string cigar);

    struct MDTagList
    {
       char operatename=0;
       int32_t size=0;
    };
    
    std::vector<MDTagList> mdtaglists;
    bool isNucleotide(char n);

    MDTagList convertMDTagList(std::string number,char operatorname);

    struct MDTagListWithSeq
    {
       char operatename=0;
       int32_t size=0;
       std::string seq;
    };

   std::vector<MDtagHelper::MDTagListWithSeq> getMDtagwithSeq(std::string *seq);
};

#endif