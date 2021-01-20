#ifndef TRAIN_H
#define TRAIN_H

#include "main.h"

#include <memory>

#include "TTree.h"

extern int Train(in_args input_args);

extern std::shared_ptr<TTree> ReadTreeFromFile(std::string tree_path, const char* tree_name);
extern std::map<std::string, int> GetTMVAMethods(std::string mymethod);

#endif