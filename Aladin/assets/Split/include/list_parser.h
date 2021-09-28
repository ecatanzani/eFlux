#ifndef LIST_PARSER_H
#define LIST_PARSER_H

#include <iostream>
#include <string>
#include <memory>
#include <fstream>
#include <sstream>

#include "TChain.h"

class parser
{
public:
    parser(const std::string input_list, const bool mc, const bool verbose);
    parser(const std::string input_list);
    ~parser(){};
    std::shared_ptr<TChain> GetEvtTree();
    const std::string GetSingleDataFile();
private:
    std::string parse_input_file(const std::string input_list);
    std::string single_data_file;
    std::shared_ptr<TChain> evtch;
};

#endif