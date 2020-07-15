#ifndef FILL_EVENT_HISTO_H
#define FILL_EVENT_HISTO_H

#include <vector>

#include "TH1D.h"
#include "TH2D.h"

#include <data_cuts.h>

extern void fill_XTRL_histo(
    const double sumRms,
    const std::vector<double> fracLayer,
    const cuts_conf data_cuts,
    const double energy,
    TH1D &h_xtrl_energy_int,
    TH2D &h_xtrl);

#endif