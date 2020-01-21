#ifndef TASK_H
#define TASK_H
#include "samplestat.h"
#include "filemanager.h"
#include <string>

class Task
{
private:
  hts_idx_t *bam_index = NULL;

  std::string taskid;
  SampleStat samplestat;
  FileManager *filepath;
  std::string target_chromosome;
  bam_hdr_t bam_header;

  struct TaskRange
  {
    uint64_t pos;
    uint64_t end;
  };

  TaskRange taskrange;

public:
  Task(SampleStat samplepath_T, FileManager *filepath, std::string target_chromosome);
  void setRange(uint64_t pos,uint64_t end);
  void setTaskName();
  void execute();
  void setBamHeader(bam_hdr_t bam_header);
  void setHtsIndex(hts_idx_t *index);
  void showTaskDone();
};

#endif