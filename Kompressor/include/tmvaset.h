#ifndef TMVASET_H
#define TMVASET_H

#include <string>
#include <memory>

#include "TChain.h"

#include "config.h"
#include "energy_config.h"

extern void createTMVAset(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    const std::string fit_tree_path,
    const bool signal,
    const std::string output_file,
    const bool _VERBOSE,
    const bool no_split,
    const bool no_split_test,
    const unsigned int threads,
    const bool _mc);

#endif