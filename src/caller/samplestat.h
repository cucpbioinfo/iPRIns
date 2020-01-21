#ifndef SAMPLESTAT_H
#define SAMPLESTAT_H
#include <htslib/sam.h>
#include <string>
#include <vector>

class SampleStat
{
private:
  std::string sample_path;

  int32_t insertsize_avg=0;
  int32_t insertsize_sd=0;
  int avg_rd=0;
  int read_length=0;
  
  //config
  int countMax=100000;

  void findReadLength();
  int32_t findMedianSampleStat(std::vector<int32_t> *insertlist,int32_t limitsize,int32_t minsize);
  int32_t findSDSampleStat(std::vector<int32_t> *insertlist,int32_t insertsize_avg,int32_t limitsize,int32_t minsize);
  std::vector<int32_t> getInsertSizeList(int64_t numberofread);

public:
  SampleStat(std::string samplepath);
  SampleStat();
  void execute();
  void setSamplePath(std::string t_sample_path);
  void setNumberOfRead(int number);
  int32_t getReadLength();
  int32_t getAverageSampleStat();
  int32_t getSDSampleStat();
};

#endif