#ifndef BGOFIDUCIAL_H
#define BGOFIDUCIAL_H

#include <memory>

#include "TFile.h"
#include "TChain.h"

extern void bgofiducial_distributions(const std::string output_path, std::shared_ptr<TChain> evtch, const bool verbose);

#endif