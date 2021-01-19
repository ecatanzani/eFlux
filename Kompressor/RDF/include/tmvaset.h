#ifndef TMVASET_H
#define TMVASET_H

#include <string>
#include <memory>

#include "TChain.h"

#include "config.h"

extern void createTMVAset(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    const std::string fit_tree_path,
    const bool signal,
    const std::string output_file,
    const bool _VERBOSE);

#endif