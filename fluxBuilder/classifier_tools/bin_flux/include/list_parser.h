#ifndef LIST_PARSER_H
#define LIST_PARSER_H

#include <string>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TChain.h"

class parser {
public:
    parser(const std::string input_list, const bool verbose);
    ~parser(){};
    std::shared_ptr<TChain> GetEvtTree();
    const std::string GetSingleDataFile();
private:
    std::string parse_input_file(const std::string input_list);
    std::shared_ptr<TChain> evtch;
};

#endif