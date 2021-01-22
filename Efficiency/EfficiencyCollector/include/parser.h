#ifndef PARSER_H
#define PARSER_H

#include <iostream>
#include <string>
#include <memory>
#include <fstream>
#include <sstream>

#include "TChain.h"

class parser
{
public:
    parser(
        const std::string input_list,
        const bool mc = false,
        const bool _VERBOSE = true);
    ~parser(){};
    std::shared_ptr<TChain> GetEvtTree();
private:
    std::string parse_input_file(const std::string input_list);
    std::string mc_tree_name = "DmpMCEvtNtup";
    std::string data_tree_name = "DmpEvtNtup";
    std::shared_ptr<TChain> evtch;
};

#endif