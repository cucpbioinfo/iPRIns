#ifndef FILEPATHMANAGER_H
#define FILEPATHMANAGER_H

#include <string>
#include <mutex>

class FileManager {
private:
    std::string sample_path;
    std::string reference_path;
    std::string output_path;
    std::string evidence_path;

    std::mutex mutexReadDepthStat;


public:
    FileManager(std::string t_sample_path, std::string t_reference_path, std::string t_output_path);

    FileManager() {};

    std::string getSamplePath();

    std::string getReferencePath();

    std::string getOutputPath();

    std::string getEvidencePath();
    std::string getVariantPath();

    std::string getResultVcfPath();

    std::string getReadDepthPath();

    std::string getAllEvidencePath();

    std::string getTempEvidencePath();

    std::string getSplitReadPath();

    std::string getReadDepthStatPath();

    void setSamplePath(std::string samplepath);

    void setReferencePath(std::string referencepath);

    void setOutputPath(std::string outputpath);


    void createDirectory(std::string yourpath);

    void initialize();

    void mutexLockReadDepthStat();
    void mutexUnlockReadDepthStat();
};

#endif