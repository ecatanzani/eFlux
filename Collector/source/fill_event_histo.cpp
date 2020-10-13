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
	auto energy_gev = energy * _GeV;

	if (lastFracLayer != -1)
	{
		auto xtrl = 0.125e-6 * pow(sumRms, 4) * lastFracLayer;
		
		// Fill XTRL histos
		h_xtrl_energy_int.Fill(xtrl, energy_w);
		h_xtrl.Fill(energy_gev, xtrl, energy_w);
		auto xtrl_bin_idx = h_xtrl.GetXaxis()->FindBin(energy_gev) -1;
		if (xtrl_bin_idx>=0 && xtrl_bin_idx<bin_xtrl.size())
			bin_xtrl[xtrl_bin_idx]->Fill(xtrl, energy_w);
	} 
}

void fill_ep_histos(
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
	const double energy_w)
{
	double _GeV = 0.001;
	auto energy_gev = energy * _GeV;

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
	
		if (energy_gev >= 20 && energy_gev < 100)
		{
			e_discrimination_20_100.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
			e_discrimination_last_20_100.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
		}
		else if (energy_gev >= 100 && energy_gev < 250)
		{
			e_discrimination_100_250.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
			e_discrimination_last_100_250.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
		}
		else if (energy_gev >= 250 && energy_gev < 500)
		{
			e_discrimination_250_500.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
			e_discrimination_last_250_500.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
		}
		else if (energy_gev >= 500 && energy_gev < 1000)
		{
			e_discrimination_500_1000.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
			e_discrimination_last_500_1000.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
		}
		else if (energy_gev >= 1000 && energy_gev < 3000)
		{
			e_discrimination_1000_3000.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
			e_discrimination_last_1000_3000.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
		}
		else
		{
			e_discrimination_3000_10000.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
			e_discrimination_last_3000_10000.Fill(
				sumRMS, 
				lastFracLayer,
				energy_w);
		}
	}
}

void fill_sumRms_cosine_histo(
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
	const double energy,
	const double energy_w)
{	
	double _GeV = 0.001;
	auto energy_gev = energy * _GeV;
	
	for (auto it=logEBins.begin(); it != logEBins.end() -1; ++it)
		if (energy_gev >= (*it) && energy_gev < (*it+1))
		{
			auto energy_idx = std::distance(logEBins.begin(), it);
			sumRms_cosine[energy_idx]->Fill(costheta, sumRMS, energy_w);
			break;
		}
	
	if (energy_gev >= 20 && energy_gev < 100)
		sumRms_cosine_20_100.Fill(costheta, sumRMS, energy_w);
	else if (energy_gev >= 100 && energy_gev < 250)
		sumRms_cosine_100_250.Fill(costheta, sumRMS, energy_w);
	else if (energy_gev >= 250 && energy_gev < 500)
		sumRms_cosine_250_500.Fill(costheta, sumRMS, energy_w);
	else if (energy_gev >= 500 && energy_gev < 1000)
		sumRms_cosine_500_1000.Fill(costheta, sumRMS, energy_w);
	else if (energy_gev >= 1000 && energy_gev < 3000)
		sumRms_cosine_1000_3000.Fill(costheta, sumRMS, energy_w);
	else
		sumRms_cosine_3000_10000.Fill(costheta, sumRMS, energy_w);
	
}