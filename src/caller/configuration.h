#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "configuration.h"
#include <stdint.h>
#include <string>

class Configuration
{
private:
    

public:
    Configuration();
    std::string mode;

    void setMode(std::string mode);
    void getMode(std::string mode);
};

#endif