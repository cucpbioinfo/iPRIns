#include <iostream>
#include <stdio.h>
#include "cli.h"

int main(int argc, char **argv)
{
    
    Cli cli(argc, argv);

    if (cli.getCommand()=="call")
    {
        return cli.callSV();
    }
    else if (cli.getCommand()=="debug")
    {
        return cli.debug();
    }
    else if (cli.getCommand()=="version")
    {
        std::cout << "VERSION:" << std::endl;
        std::cout << "\t1.1.0" << std::endl;
        return 0;
    }
    else if (cli.getCommand()=="")
    {
        cli.ShowHelp();
        return 0;
    }
    else
    {
        std::cout << cli.getCommand() <<" : command not found" << std::endl;
        std::cout << std::endl;
        cli.ShowHelp();
        return 1;
    }

    return 0;
}