#include "filemanager.h"
#include <iostream>
#include <sys/stat.h>
#include <sstream>

FileManager::FileManager(std::string t_sample_path, std::string t_reference_path, std::string t_output_path)
{
    sample_path = t_sample_path;
    reference_path = t_reference_path;
    output_path = t_output_path;

    // std::string s;
    // s.append(*output_path);
    // s.append("/analysis/evidence");
    evidence_path = output_path + "/analysis/evidence";
}

std::string FileManager::getOutputPath()
{
    return output_path;
}

std::string FileManager::getAllEvidencePath()
{
    return getOutputPath()+"/analysis/" + "allevidence.txt";
}

std::string FileManager::getTempEvidencePath()
{
    return getOutputPath()+"/analysis/" + "evidence";
}

std::string FileManager::getSplitReadPath() {
    return getOutputPath()+"/analysis/" + "splitread";
}


std::string FileManager::getVariantPath()
{
    return getOutputPath()+"/analysis/" + "variant/";
}

std::string FileManager::getResultVcfPath() {
    return getOutputPath()+"/result.vcf";
}

std::string FileManager::getEvidencePath()
{
    return getOutputPath()+"/analysis/" + "evidence";
}

std::string FileManager::getSamplePath()
{
    return sample_path;
}

std::string FileManager::getReferencePath()
{
    return reference_path;
}

void FileManager::setSamplePath(std::string T_sample_path)
{
    sample_path = T_sample_path;
}

void FileManager::setReferencePath(std::string T_reference_path)
{
    reference_path = T_reference_path;
}


void FileManager::setOutputPath(std::string T_output_path)
{
    output_path = T_output_path;
}

void FileManager::initialize()
{
    createDirectory(output_path);
    createDirectory(output_path + "/analysis");
    createDirectory(output_path + "/analysis/evidence");
    createDirectory(output_path + "/analysis/readdepth");
    createDirectory(output_path + "/analysis/readdepthstat");
    createDirectory(output_path + "/analysis/variant");
    createDirectory(output_path + "/analysis/splitread");
}

std::string FileManager::getReadDepthStatPath()
{
    return getOutputPath()+"/analysis/readdepthstat/readdepthstat.txt";
}

std::string FileManager::getReadDepthPath()
{
    return output_path + "/analysis/readdepth";
}

void FileManager::createDirectory(std::string yourpath)
{
    if (mkdir(yourpath.c_str(), 0777) == -1)
    {
        // std::cerr << "Error :  " << strerror(errno) << std::endl;
        // exit(0);
    }
    // std::cerr << "create directory : " << yourpath << " = " << strerror(errno) << std::endl;
}

void FileManager::mutexLockReadDepthStat()
{
    mutexReadDepthStat.lock();
}

void FileManager::mutexUnlockReadDepthStat()
{
    mutexReadDepthStat.unlock();
}
