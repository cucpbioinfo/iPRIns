#ifndef SPECIFYINGEVIDENCEINSERTION_H
#define SPECIFYINGEVIDENCEINSERTION_H
#include "specifyingevidence.h"
#include "evidence.h"

class SpecifyingEvidenceInsertion : public SpecifyingEvidence
{
private:
  void checkRange();

  uint32_t currentPos = 0;
  uint32_t currentMPos = 0;
  void proveEvidence(int index);
  void checkProveEvidence();
  bool filterEvidence(Evidence *evidence);
  void calculateVCF(Evidence *evidence);
protected:

public:
  SpecifyingEvidenceInsertion();
  void updateRead();
  void done();
  bool incrementSVFreq(int32_t overlappedpos, int32_t pos, int32_t mpos);
};

#endif