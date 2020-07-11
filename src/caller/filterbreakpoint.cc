#include "filterbreakpoint.h"
#include <dirent.h>

FilterBreakpoint::FilterBreakpoint() {
}

int32_t FilterBreakpoint::roundNumber(int32_t number, int32_t round) {
    if (number % round == 0) {
        return number;
    }

    int32_t tempNumber = 0;
    if (number % round >= round / 2) {
        int32_t numberrounddiff = round - (number % round);
        tempNumber = numberrounddiff + number;
    } else {
        int32_t numberrounddiff = (number % round);
        tempNumber = number - numberrounddiff;
    }

    return tempNumber;
}

int32_t FilterBreakpoint::nextNumber(int32_t number, int32_t round) {
    return roundNumber(number, round) + round;
}

int32_t FilterBreakpoint::previousNumber(int32_t number, int32_t round) {
    return roundNumber(number, round) - round;
}

void FilterBreakpoint::execute() {
    std::vector <std::string> vcffilelist = getPathVCFFiles();

    for (auto n : vcffilelist) {
        auto variantlist = getEvidenceByFilepath(n);

        if (variantlist.size() == 0) {
            continue;
        }

        std::vector <Evidence> result;
        if (variantlist.at(0).getSVType() == "INS") {

            // result = getResultWithOutOverlapped(&variantlist, &variantlist);
            result = getRefineResultInsertion(&variantlist);
            result = getResultRemoveOverlapped(&result, &result);

            // result = variantlist;
        } else {
            result = variantlist;
        }

        writeFile(&result);
    }
}

void FilterBreakpoint::setSampleStat(SampleStat *samplestat) {
    FilterBreakpoint::samplestat = samplestat;
}

std::vector <Evidence> FilterBreakpoint::getRefineResultInsertion(std::vector <Evidence> *master) {

    std::vector <Evidence> cache;
    for (auto n : *master) {

        if ((n.getMark() == "MATEUNMAPPED")) {

            if (n.getFrequency() > readDepthStat.getReadDepthByChr(n.getChr()) * 2) {
                continue;
            }

            if (n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 10, 100, 4)) {
                continue;
            }

            if (n.getNumberOfRP() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 10, 100, 4)) {
                continue;
            }

            if (n.LNGMATCH <= getDivider(samplestat->getReadLength(), 30, 100, 15)) {
                continue;
            }

            if (n.getMaxMapQ() < 20) {
                continue;
            }

            if (n.getMaxRPMapQ() < 40) {
                continue;
            }

            if (n.getMaxMapQ() >= 60 || n.getMaxRPMapQ() >= 60) {
            } else {
                continue;
            }

        } else if ((n.getMark() == "SINS")) {
            if (n.getSvLength() < 50) {
                continue;
            }

        } else if ((n.getMark() == "SR")) {
            if (n.getNumberOfRP() <= 1) {
                continue;
            }

        } else if ((n.getMark() == "UNMERGE")) {
            if (n.getFrequency() > readDepthStat.getReadDepthByChr(n.getChr()) * 3) {
                continue;
            }

            if (n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 2, 100, 3)) {
                continue;
            }
//            continue;
        } else {
//            continue;

            if (n.getFrequency() > readDepthStat.getReadDepthByChr(n.getChr()) * 3) {
                continue;
            }

//            if (n.getMaxMapQ() < 10)
//            {
//                continue;
//            }
//
//            if (n.getFrequency() <= 2)
//            {
//                continue;
//            }
//
//            if (n.getNumberOfRP() <= 3)
//            {
//                continue;
//            }

//            if (n.getFrequency() <= 4)
//            {
//                if (n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 2))
//                {
//                    continue;
//                }
//            }
//
//            if (n.getNumberOfRP() <= 5)
//            {
//                if (n.getNumberOfRP() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 7, 100, 3))
//                {
//                    continue;
//                }
//            }

            if (n.getMaxMapQ() < 60) {
                if (n.LNGMATCH <= getDivider(samplestat->getReadLength(), 15, 100, 15)) {
                    continue;
                }
            } else {
                if (n.LNGMATCH <= getDivider(samplestat->getReadLength(), 20, 100, 15)) {
                    continue;
                }
            }

            if (n.getMaxMapQ() >= 60 || n.getMaxRPMapQ() >= 60) {
            } else {
                continue;
            }

            // continue;
        }

        cache.push_back(n);
    }

    return cache;
}


void FilterBreakpoint::writeFile(std::vector <Evidence> *master) {
    std::ofstream myfile;
    myfile.open(filemanager->getOutputPath() + "/result.vcf", std::ios_base::app);
    for (auto n : *master) {
        n.setID("iPRIns" + std::to_string(vcfIdNumber));
        myfile << n.getResultVcfFormatString() << std::endl;
        vcfIdNumber++;
    }

    myfile.close();
}

std::vector <Evidence>
FilterBreakpoint::getResultWithOutOverlapped(std::vector <Evidence> *master, std::vector <Evidence> *slave) {
    std::vector <Evidence> cache;
    for (auto n : *master) {
        bool found = false;
        for (auto m : *slave) {
            if (n.getPos() == m.getPos() && n.getEnd() == m.getEnd()) {
                continue;
            }

            if (m.getPos() - 100 <= n.getPos() && n.getPos() <= m.getPos() + 100) {
                found = true;
                break;
            }
        }

        if (!found) {
            cache.push_back(n);
        }
    }

    return cache;
}

std::vector <Evidence>
FilterBreakpoint::getResultRemoveOverlapped(std::vector <Evidence> *master, std::vector <Evidence> *slave) {
    std::vector <Evidence> cache;
    for (auto n : *master) {
        bool found = false;
        for (auto m : cache) {
            if (n.getPos() == m.getPos() && n.getEnd() == m.getEnd()) {
                found = true;
                break;
            }

            if (checkBetween(n.getPos(), m.getPos(), 10) && checkBetween(n.getEnd(), m.getEnd(), 10)) {
                // std::cout << n.getChr() << " " << n.getPos() << " " << m.getEnd() << std::endl;
                found = true;
                break;
            }
        }

        if (!found) {
            cache.push_back(n);
        }
    }

    return cache;
}

bool FilterBreakpoint::checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped) {
    if (targetPos - overlapped > pos) {
        return false;
    }

    if (targetPos + overlapped < pos) {
        return false;
    }

    return true;
}

std::vector <Evidence> FilterBreakpoint::getEvidenceByFilepath(std::string filepaht) {
    std::vector <Evidence> cache;
    std::string line;
    std::ifstream myfile(filepaht);
    if (myfile.is_open()) {
        while (getline(myfile, line)) {
            Evidence e;
            e.setEvidenceByString(line);
            if (e.getPos() != 0) {
                cache.push_back(e);
            }

        }
        myfile.close();
    }

    return cache;
}

void FilterBreakpoint::setFileManager(FileManager *filemanager) {
    FilterBreakpoint::filemanager = filemanager;
    readDepthStat.setFilePath(filemanager);
    readDepthStat.execute();
}

std::vector <std::string> FilterBreakpoint::getPathVCFFiles() {
    std::vector <std::string> evidenceFilePathLists;
    DIR *d;
    struct dirent *dir;

    d = opendir(filemanager->getVariantPath().c_str());
    if (d) {
        while (dir = readdir(d)) {
            if (std::string(dir->d_name).size() < 4) {
                continue;
            }

            if (std::string(dir->d_name).substr(std::string(dir->d_name).size() - 4) == ".vcf") {
                std::string tempPath = filemanager->getVariantPath() + std::string(dir->d_name);
                // std::cout << "tempPath : " << tempPath << std::endl;
                evidenceFilePathLists.push_back(tempPath);
            }
        }
        closedir(d);
    }

    return evidenceFilePathLists;
}

int FilterBreakpoint::getDivider(int value, int top, int down, int minimum) {
    auto returnvalue = (int) (float(value) * (float(top) / float(down)));

    if (returnvalue > minimum) {
        return returnvalue;
    }

    return minimum;
}