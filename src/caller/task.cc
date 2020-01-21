
#include "task.h"
#include <map>
#include <iostream>
#include "evidencefinder.h"

Task::Task(SampleStat samplestat_T, FileManager *filepath_T, std::string target_chromosome_T)
{
    samplestat = samplestat_T;
    filepath = filepath_T;
    target_chromosome = target_chromosome_T;
}

void Task::setRange(uint64_t t_pos, uint64_t t_end)
{
    taskrange.pos = t_pos;
    taskrange.end = t_end;
}

void Task::setBamHeader(bam_hdr_t t_bam_header)
{
    bam_header = t_bam_header;
}

void Task::setHtsIndex(hts_idx_t *index)
{
    bam_index = index;
}

void Task::execute()
{
    EvidenceFinder evidencefinder(&samplestat, filepath, &target_chromosome);
    evidencefinder.setBamHeader(&bam_header);
    evidencefinder.setHtsIndex(bam_index);
//     std::string compare = "chr1";
//     if (target_chromosome == compare)
//     {
        evidencefinder.execute();
//     }
}

void Task::showTaskDone()
{
    std::cout << "✓ " << target_chromosome << std::endl;
}
