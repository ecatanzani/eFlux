#ifndef FIT_H
#define FIT_H

#include <memory>

#include "config.h"
#include "energy_config.h"
#include "lambda_config.h"

#include "TChain.h"

extern void fit(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    std::shared_ptr<lambda_config> _lambda_config,
    const double _entries,
    unsigned int focus_energybin,
    const std::string outputPath,
    const std::string regularize_tree_path,
    const bool verbose,
    const unsigned int threads,
    const bool _mc);

#endif