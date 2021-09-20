#ifndef KOMPRESS_H
#define KOMPRESS_H

#include <string>
#include <memory>

#include "config.h"
#include "energy_config.h"

#include "TChain.h" 

extern void kompress(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    const std::string fit_tree_path,
    const double _entries,
    const std::string outputPath,
    const bool verbose,
    const unsigned int threads,
    const bool _mc);

#endif