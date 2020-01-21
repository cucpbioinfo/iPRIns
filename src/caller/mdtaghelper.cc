#include "mdtaghelper.h"
#include <iostream>

MDtagHelper::MDtagHelper(std::string mdtext,std::string cigar)
{

    std::string mergeNumber;
    char operatorname;
    int sizeOfMDtext = mdtext.size();
    int count=0;

    for (auto n : mdtext)
    {
        count++;
        if (isdigit(n))
        {
            mergeNumber += n;
            if (count==sizeOfMDtext)
            {
                mdtaglists.push_back(convertMDTagList(mergeNumber,'M'));
            }
        }
        else
        {
            if (isNucleotide(n))
            {
                if (mergeNumber != "")
                {
                    mdtaglists.push_back(convertMDTagList(mergeNumber,'M'));
                    operatorname = '\0';
                    mergeNumber = "";
                }
                
                    mdtaglists.push_back(convertMDTagList("1", 'X'));
                    operatorname = '\0';
                    mergeNumber = "";
            }
            else if (n=='^')
            {
                mdtaglists.push_back(convertMDTagList("1", 'D'));
                operatorname = '\0';
                mergeNumber = "";
            }
            else
            {
                exit(0);
            }
        }
    }

}

MDtagHelper::MDTagList MDtagHelper::convertMDTagList(std::string number, char operatorname)
{
    MDTagList md;
    md.size = std::stoi(number);
    md.operatename = operatorname;
    return md;
}

bool MDtagHelper::isNucleotide(char n)
{
    switch (n)
    {
    case 'A':
        return true;
    case 'C':
        return true;
    case 'T':
        return true;
    case 'G':
        return true;
    default: //Optional
        return false;
    }

    return false;
}

std::vector<MDtagHelper::MDTagListWithSeq> MDtagHelper::getMDtagwithSeq(std::string *seq)
{
    std::vector<MDTagListWithSeq> result;
    std::string sumText;
    
    int currentPrevious;
    std::cout << "----------------------------" << mdtaglists.size() << std::endl;
    int currentSize = (*seq).size();
    int currentPos = currentSize;
    for (int i=mdtaglists.size()-1;i>-1;i--)
    {
            currentPos = currentPos-mdtaglists[i].size;

            // std::string textCurrent = (*seq).substr(currentSize-mdtaglists[i].size,mdtaglists[i].size);
            // sumText.insert(0,textCurrent);
            // std::cout << sumText << std::endl;

                // std::cout << mdtaglists[i].operatename << " / " << mdtaglists[i].size << std::endl;
                // std::cout << sumText << std::endl;
                // std::string stringInsert = "";
                
                // sumText.insert(0,);
                // str.insert(0,str2);
                // (*seq).substr(1,(*seq).size());

    }
   
    return result;
}