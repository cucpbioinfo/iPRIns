#include "variantresultfilter.h"

VariantResultFilter::VariantResultFilter()
{
}

bool VariantResultFilter::passFilterSV(Evidence *variantresult)
{
    if (!variantresult->isQuailtyPass())
    {
        return false;
    }


    if (variantresult->getChr() == "")
    {
        return false;
    }

    if (variantresult->getPos() == 0)
    {
        return false;
    }

    return true;
}

bool VariantResultFilter::passFilterInsertion(Evidence *variantresult)
{
    if (variantresult->getFrequency() <= 2)
    {
        return false;
    }

    return true;
}
