#ifndef FILL_EVENT_HISTO_H
#define FILL_EVENT_HISTO_H

#include <vector>
#include <memory>

#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"

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
	TH2D &e_discrimination_20_100,
	TH2D &e_discrimination_100_250,
	TH2D &e_discrimination_250_500,
	TH2D &e_discrimination_500_1000,
	TH2D &e_discrimination_1000_3000,
	TH2D &e_discrimination_3000_10000,
	TH2D &e_discrimination_last,
	TH2D &e_discrimination_last_20_100,
	TH2D &e_discrimination_last_100_250,
	TH2D &e_discrimination_last_250_500,
	TH2D &e_discrimination_last_500_1000,
	TH2D &e_discrimination_last_1000_3000,
	TH2D &e_discrimination_last_3000_10000,
	const double energy,
	const double energy_w = 1);

extern void fill_sumRms_cosine_histo(
	const double sumRMS,
	const double costheta,
	const std::vector<float> logEBins,
	std::vector<std::shared_ptr<TH2D>> &sumRms_cosine,
	TH2D &sumRms_cosine_20_100,
	TH2D &sumRms_cosine_100_250,
	TH2D &sumRms_cosine_250_500,
	TH2D &sumRms_cosine_500_1000,
	TH2D &sumRms_cosine_1000_3000,
	TH2D &sumRms_cosine_3000_10000,
	const double energy_raw,
	const double energy_w = 1);

#endif