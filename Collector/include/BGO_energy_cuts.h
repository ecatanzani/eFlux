#ifndef BGO_ENERGY_CUTS_H
#define BGO_ENERGY_CUTS_H

#include <memory>
#include <vector>
#include <algorithm>

#include "DmpEvtBgoRec.h"

#include "TH1D.h"

#include "data_cuts.h"
#include "DAMPE_geo_structure.h"

#if 0
extern void evaluateEnergyRatio(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const double bgoTotalE,
    TH1D &h_layer_max_energy_ratio,
    std::vector<TH1D> &h_layer_energy_ratio,
    const double energy_w = 1);
#endif

extern void evaluateEnergyRatio(
    const std::vector<double> fracLayer,
    TH1D &h_layer_max_energy_ratio,
    std::vector<TH1D> &h_layer_energy_ratio,
    const double energy_w = 1);

#endif