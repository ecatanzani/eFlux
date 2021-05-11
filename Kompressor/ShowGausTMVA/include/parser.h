#ifndef PARSER_H
#define PARSER_H

#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TChain.h"

extern std::string parse_input_file(const std::string input_list);
extern std::shared_ptr<TChain> getchain(
    const std::string filelist, 
    const bool mc = false,
    const bool verbose = true);

#endif

