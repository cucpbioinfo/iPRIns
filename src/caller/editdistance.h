#ifndef EEITDISTANCE_H
#define EEITDISTANCE_H
#include <string>

class EditDistance
{
private:

public:
  EditDistance();
  int Compare(std::string *str1,std::string *str2);
};

#endif