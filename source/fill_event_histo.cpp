#include "fill_event_histo.h"

void fill_XTRL_histo(
	const double sumRms,
	const double lastFracLayer,
	const double energy,
	std::vector<std::shared_ptr<TH1D>> &bin_xtrl,
	TH1D &h_xtrl_energy_int,
	TH2D &h_xtrl)
{
	double _GeV = 0.001;

	if (lastFracLayer != -1)
	{
		auto xtrl = 0.125e-6 * pow(sumRms, 4) * lastFracLayer;
		
		// Fill XTRL histos
		h_xtrl_energy_int.Fill(xtrl);
		h_xtrl.Fill(energy * _GeV, xtrl);
		bin_xtrl[(h_xtrl.GetXaxis()->FindBin(energy * _GeV)) - 1]->Fill(xtrl);
	} 
}

void fill_ep_histos(
	const double sumRMS,
	const double lastFracLayer,
	const double frac_layer_13,
	TH2D &e_discrimination,
	TH2D &e_discrimination_last)
{
	if (lastFracLayer != -1 && frac_layer_13 != -1)
	{
		e_discrimination.Fill(
			sumRMS, 
			lastFracLayer);
		e_discrimination_last.Fill(
			sumRMS, 
			frac_layer_13);
	}
}