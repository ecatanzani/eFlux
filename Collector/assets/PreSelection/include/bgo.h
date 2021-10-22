#ifndef BGO_H
#define BGO_H

#include <memory>

#include "main.h"
#include "config.h"

#include "TFile.h"
#include "TChain.h"

extern void bgofiducial_distributions(in_pars input_pars, std::shared_ptr<TChain> evtch);
extern void bgofiducial_distributions_mc(in_pars input_pars, std::shared_ptr<TChain> evtch);
extern void bgofiducial_distributions_data(in_pars input_pars, std::shared_ptr<TChain> evtch);

#endif