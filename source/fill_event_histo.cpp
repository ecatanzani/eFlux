#include "fill_event_histo.h"

void fill_XTRL_histo(
	const double sumRms,
	const double lastFracLayer,
	const double energy,
	std::vector<std::shared_ptr<TH1D>> &bin_xtrl,
	TH1D &h_xtrl_energy_int,
	TH2D &h_xtrl,
	const double energy_w)
{
	double _GeV = 0.001;

	if (lastFracLayer != -1)
	{
		auto xtrl = 0.125e-6 * pow(sumRms, 4) * lastFracLayer;
		
		// Fill XTRL histos
		h_xtrl_energy_int.Fill(xtrl, energy_w);
		h_xtrl.Fill(energy * _GeV, xtrl, energy_w);
		bin_xtrl[(h_xtrl.GetXaxis()->FindBin(energy * _GeV)) - 1]->Fill(xtrl, energy_w);
	} 
}

void fill_ep_histos(
	const double sumRMS,
	const double lastFracLayer,
	const double frac_layer_13,
	TH2D &e_discrimination,
	TH2D &e_discrimination_last,
	const double energy_w)
{
	if (lastFracLayer != -1 && frac_layer_13 != -1)
	{
		e_discrimination.Fill(
			sumRMS, 
			lastFracLayer,
			energy_w);
		e_discrimination_last.Fill(
			sumRMS, 
			frac_layer_13,
			energy_w);
	}
}

void fill_sumRms_cosine_histo(
	const double sumRMS,
	const double costheta,
	const double energy_corr,
	const std::vector<float> logEBins,
	std::vector<std::shared_ptr<TH2D>> &sumRms_cosine,
	TH2D &sumRms_cosine_20_100,
	TH2D &sumRms_cosine_100_250,
	TH2D &sumRms_cosine_250_500,
	TH2D &sumRms_cosine_500_1000,
	TH2D &sumRms_cosine_1000_3000,
	TH2D &sumRms_cosine_3000_10000)
{
	for (auto it=logEBins.begin(); it != logEBins.end() -1; ++it)
		if (energy_corr >= (*it) && energy_corr < (*it+1))
		{
			auto energy_idx = std::distance(logEBins.begin(), it);
			sumRms_cosine[energy_idx]->Fill(costheta, sumRMS);
			break;
		}
	
	if (energy_corr >= 20 && energy_corr < 100)
		sumRms_cosine_20_100.Fill(costheta, sumRMS);
	else if (energy_corr >= 100 && energy_corr < 250)
		sumRms_cosine_100_250.Fill(costheta, sumRMS);
	else if (energy_corr >= 250 && energy_corr < 500)
		sumRms_cosine_250_500.Fill(costheta, sumRMS);
	else if (energy_corr >= 500 && energy_corr < 1000)
		sumRms_cosine_500_1000.Fill(costheta, sumRMS);
	else if (energy_corr >= 1000 && energy_corr < 3000)
		sumRms_cosine_1000_3000.Fill(costheta, sumRMS);
	else
		sumRms_cosine_3000_10000.Fill(costheta, sumRMS);
	
}