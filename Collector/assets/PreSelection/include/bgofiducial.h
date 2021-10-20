#ifndef BGOFIDUCIAL_H
#define BGOFIDUCIAL_H

#include <memory>

#include "config.h"

#include "TFile.h"
#include "TChain.h"

extern void bgofiducial_distributions(
    const std::string output_path, 
    const std::string logs_dir, 
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> evt_config,
    const bool verbose);

#endif