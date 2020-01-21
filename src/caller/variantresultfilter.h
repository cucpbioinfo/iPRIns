#ifndef VARIANTRESULTFILTER_H
#define VARIANTRESULTFILTER_H
#include "evidence.h"
class VariantResultFilter
{
private:
    
    
public:
    VariantResultFilter();
    bool passFilterSV(Evidence *variantresult);
    bool passFilterInsertion(Evidence *variantresult);

};

#endif