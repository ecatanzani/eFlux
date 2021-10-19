#ifndef BGOFIDUCIAL_H
#define BGOFIDUCIAL_H

#include <memory>

#include "TFile.h"
#include "TChain.h"

extern void bgofiducial_distributions(TFile* outfile, std::shared_ptr<TChain> evtch, const bool verbose);

#endif