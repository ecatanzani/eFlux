#ifndef FILL_EVENT_HISTO_H
#define FILL_EVENT_HISTO_H

#include <vector>
#include <memory>

#include "TH1D.h"
#include "TH2D.h"

#include <data_cuts.h>

extern void fill_XTRL_histo(
	const double sumRms,
	const double lastFracLayer,
	const double energy,
	std::vector<std::shared_ptr<TH1D>> &bin_xtrl,
	TH1D &h_xtrl_energy_int,
	TH2D &h_xtrl,
	const double energy_w = 1);

extern void fill_ep_histos(
	const double sumRMS,
	const double lastFracLayer,
	const double frac_layer_13,
	TH2D &e_discrimination,
	TH2D &e_discrimination_last,
	const double energy_w = 1);

#endif