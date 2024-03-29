#ifndef READER_H
#define READER_H

#include <string>
#include <memory>

#include "TChain.h"

extern std::shared_ptr<TChain> BuildChain(const std::string input_list, const bool verbose);
extern std::shared_ptr<TChain> ReadTreeFromFile(const std::string input_list, const char *tree_name, const bool verbose=false);
extern std::string parse_input_file(std::string input_list);

#endif