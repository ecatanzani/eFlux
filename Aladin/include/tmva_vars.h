#ifndef TMVA_VARS_H
#define TMVA_VARS_H

#include "config.h"
#include "energy_config.h"

#include "TChain.h"

extern void tmva_vars(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<energy_config> _energy_config,
    const double _entries,
    const std::string lambda_tree,
    const std::string regularize_tree,
    const std::string outputPath,
    const bool verbose,
    const unsigned int threads,
    const bool _mc);

#endif