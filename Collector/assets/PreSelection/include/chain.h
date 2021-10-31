#ifndef CHAIN_H
#define CHAIN_H

#include "TChain.h"

#include <memory>
#include <string>

extern std::shared_ptr<TChain> GetFileChain(const std::string input, const bool verbose, const bool mc);

#endif