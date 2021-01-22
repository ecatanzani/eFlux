#ifndef READER_H
#define READER_H

#include <memory>
#include <string>
#include <memory>

#include "TChain.h" 

#include "main.h"
#include "config.h"
#include "energy_config.h"

extern void reader(in_args input_args);

extern void mc_reader(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    const std::string fit_tree_path,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE,
    const unsigned int threads);

#endif
