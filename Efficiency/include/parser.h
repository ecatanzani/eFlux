#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TChain.h"

class parser
{
public:
    parser(const std::string input_list, const bool verbose);
    ~parser(){};
    std::shared_ptr<TChain> GetEvtTree();
    const std::string GetSingleDataFile();
private:
    std::string parse_input_file(const std::string input_list);
    std::string tree_name = "DmpMCEvtNtup";
    std::shared_ptr<TChain> evtch;
};

#endif