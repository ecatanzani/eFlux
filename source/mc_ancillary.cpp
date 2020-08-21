#include "mc_ancillary.h"

void compute_proton_background(
	const double sumRms,
	const double lastFracLayer,
	const double energy,
	const cuts_conf data_cuts,
	TH1D &h_background_under_xtrl_cut,
	TH1D &h_background_over_xtrl_cut,
	const double energy_w)
{
	double _GeV = 0.001;

	if (lastFracLayer != -1)
	{
		auto xtrl = 0.125e-6 * pow(sumRms, 4) * lastFracLayer;	
		if (xtrl < data_cuts.xtrl)
			h_background_under_xtrl_cut.Fill(energy * _GeV, energy_w);
		if (xtrl > 20 && xtrl < 100)
			h_background_over_xtrl_cut.Fill(energy * _GeV, energy_w);
	}
}