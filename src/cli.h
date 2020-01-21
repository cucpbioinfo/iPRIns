#ifndef CLI_H
#define CLI_H
#include <string>
#include <vector>

class Cli
{
  private:
    int argc;
    char **argv;

    std::string name;
    std::string description;
    std::string version;
    std::vector<std::string> args_lists;

  public:
    Cli(int argc, char **argv);
    std::string getCommand();
    int callSV();
    void showHelpCallSV();
    int debug();
    void ShowHelp();

};

#endif