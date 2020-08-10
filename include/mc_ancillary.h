#ifndef MC_ANCILLARY_H
#define MC_ANCILLARY_H

#include <vector>

#include "data_cuts.h"

extern void compute_proton_background(
	const double sumRms,
	const double lastFracLayer,
	const double energy,
	const cuts_conf data_cuts,
	TH1D &h_background_under_xtrl_cut,
	TH1D &h_background_over_xtrl_cut,
	const double energy_w = 1);

#endif