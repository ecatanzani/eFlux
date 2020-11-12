#include "histos.h"
#include "binning.h"
#include "DAMPE_geo_structure.h"

#include "TMath.h"
#include "TDirectory.h"

histos::histos(std::vector<float> energy_bins)
{
	logEBins = energy_bins;
	h_BGOrec_layer_max_energy_ratio.resize(logEBins.size() - 1);
	h_BGOrec_layer_energy_ratio.resize(logEBins.size() - 1);
	h_BGOrec_layer_rms.resize(logEBins.size() - 1);
	h_BGOrec_sumrms.resize(logEBins.size() - 1);
	h_BGOrec_sumrms_weighted.resize(logEBins.size() - 1);
	h_BGOrec_sumrms_cosine.resize(logEBins.size() - 1);
	h_BGOrec_fraclast.resize(logEBins.size() - 1);
	h_BGOrec_fraclast_cosine.resize(logEBins.size() - 1);
	h_BGOrec_frac13.resize(logEBins.size() - 1);
	h_BGOrec_last_layer.resize(logEBins.size() - 1);
	h_BGOrec_hits.resize(logEBins.size() - 1);
	h_BGOrec_energy_frac_1R.resize(logEBins.size() - 1);
	h_BGOrec_energy_frac_2R.resize(logEBins.size() - 1);
	h_BGOrec_energy_frac_3R.resize(logEBins.size() - 1);
	h_BGOrec_energy_frac_5R.resize(logEBins.size() - 1);
	h_BGOrec_slopeX.resize(logEBins.size() - 1);
	h_BGOrec_slopeY.resize(logEBins.size() - 1);
	h_BGOrec_interceptX.resize(logEBins.size() - 1);
	h_BGOrec_interceptY.resize(logEBins.size() - 1);
	h_BGOrec_topMap.resize(logEBins.size() - 1);
	h_BGOrec_bottomMap.resize(logEBins.size() - 1);
	sumRms_cosine.resize(logEBins.size() - 1);
	e_discrimination_last.resize(logEBins.size() - 1);
	e_discrimination.resize(logEBins.size() - 1);
	h_xtrl_bin.resize(logEBins.size() - 1);
	h_BGOrec_shower_profile.resize(logEBins.size() - 1);
	h_BGOrec_shower_profile_cosine_upto_09.resize(logEBins.size() - 1);
	h_BGOrec_shower_profile_cone_from_09.resize(logEBins.size() - 1);

	sumRms_bins = createLogBinning(10, 2e+3, 100);
	cosine_bins = createLinearBinning(0, 1, 10);
	flast_binning = createLogBinning(1e-5, 2e-1, 100);
	xtrl_bins = createLinearBinning(0, 150, 100);

	init_trigger_histos();
	init_preselection_histos();
	init_geometric_histos();
	init_BGOfiducial_histos();
	init_BGO_histos();
	init_xtrl_histos();
	init_ep_histos();
	init_STK_histos();
	init_psd_charge_histos();
	init_stk_charge_histos();
	init_nud_histos();
}

const int histos::GetEnergyBin(const double energy)
{
	int idx = h_trigger->FindBin(energy) - 1;
	if (idx < 0 || (unsigned int)idx > (logEBins.size() - 2))
		return -999;
	else
		return idx;
}

void histos::init_trigger_histos()
{
	h_trigger = std::make_unique<TH1D>(
		"h_trigger",
		"Energy Distribution of the triggered particles; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
}

void histos::init_preselection_histos()
{
	h_geometric_cut = std::make_unique<TH1D>(
		"h_geometric_cut",
		"Energy Distribution - geometric (trigger selection) cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_maxElayer_cut = std::make_unique<TH1D>(
		"h_maxElayer_cut",
		"Energy Distribution - maxElayer cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_maxBarLayer_cut = std::make_unique<TH1D>(
		"h_maxBarLayer_cut",
		"Energy Distribution - maxBarLayer cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOTrackContainment_cut = std::make_unique<TH1D>(
		"h_BGOTrackContainment_cut", "Energy Distribution - BGOTrackContainment cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGO_fiducial_cut = std::make_unique<TH1D>(
		"h_BGO_fiducial_cut", "Energy Distibution - BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_all_cut = std::make_unique<TH1D>(
		"h_all_cut", "Energy Distribution - All cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
}

void histos::init_geometric_histos()
{
	h_geometric_maxElayer_cut = std::make_unique<TH1D>(
		"h_geometric_maxElayer_cut",
		"Energy Distribution - maxElayer + geometric cut; Real Energy (GeV); counts",
		logEBins.size() - 1, &(logEBins[0]));
	h_geometric_maxBarLayer_cut = std::make_unique<TH1D>(
		"h_geometric_maxBarLayer_cut",
		"Energy Distribution - maxBarLayer + geometric cut; Real Energy (GeV); counts",
		logEBins.size() - 1, &(logEBins[0]));
	h_geometric_BGOTrackContainment_cut = std::make_unique<TH1D>(
		"h_geometric_BGOTrackContainment_cut",
		"Energy Distribution - BGOTrackContainment + geometric cut; Real Energy (GeV); counts",
		logEBins.size() - 1, &(logEBins[0]));
	h_geometric_BGO_fiducial_cut = std::make_unique<TH1D>(
		"h_geometric_BGO_fiducial_cut",
		"Energy Distibution - BGO fiducial + geometric cut; Real Energy (GeV); counts",
		logEBins.size() - 1, &(logEBins[0]));
	h_geometric_all_cut = std::make_unique<TH1D>(
		"h_geometric_all_cut",
		"Energy Distribution - All + geometric cut; Real Energy (GeV); counts",
		logEBins.size() - 1, &(logEBins[0]));
}

void histos::init_BGOfiducial_histos()
{
	h_BGOfiducial_nBarLayer13_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_nBarLayer13_cut",
		"Energy Distribution - nBarLayer13 + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_maxRms_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_maxRms_cut",
		"Energy Distribution - maxRms  + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_track_selection_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_track_selection_cut",
		"Energy Distribution - track selection + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_psd_stk_match_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_psd_stk_match_cut",
		"Energy Distribution - psd-stk match + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_psd_charge_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_psd_charge_cut",
		"Energy Distribution - psd charge + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_stk_charge_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_stk_charge_cut",
		"Energy Distribution - stk charge + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_all_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_all_cut",
		"Energy Distribution - All + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
}

void histos::init_BGO_histos()
{
	h_BGOrec_energy = std::make_unique<TH1D>(
		"h_BGOrec_energy",
		"BGO Energy: Raw Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOrec_corr_energy = std::make_unique<TH1D>(
		"h_BGOrec_corr_energy",
		"BGO Corrected Energy: Raw Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOrec_layer_energy_diff = std::make_unique<TH1D>(
		"h_BGOrec_layer_energy_diff",
		"Sum Energy Layer - Raw Energy",
		10, 0, 1);
	h_BGOrec_costheta = std::make_unique<TH1D>(
		"h_BGOrec_costheta",
		"BGO Reco costheta",
		10, 0, 1);
	h_BGOrec_costheta_presel = std::make_unique<TH1D>(
		"h_BGOrec_costheta_presel",
		"BGO Reco costheta - After Preselection Cuts",
		10, 0, 1);

	for (unsigned int bidx = 0; bidx < (logEBins.size() - 1); ++bidx)
	{
		h_BGOrec_layer_max_energy_ratio[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_layer_max_energy_ratio_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("Max Layer Energy Ratio - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
			10, 0, 1);
		h_BGOrec_shower_profile[bidx] = std::make_unique<TH2D>(
			(std::string("h_BGOrec_shower_profile_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("Shower Profile - Reco Energy bin") + std::to_string(bidx + 1)).c_str(),
			14, 0, 13,
			10, 0, 1);
		h_BGOrec_shower_profile_cosine_upto_09[bidx] = std::make_unique<TH2D>(
			(std::string("h_BGOrec_shower_profile_cosine_upto_09_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("Shower Profile (cosine < 0.9) - Reco Energy bin") + std::to_string(bidx + 1)).c_str(),
			14, 0, 13,
			10, 0, 1);
		h_BGOrec_shower_profile_cone_from_09[bidx] = std::make_unique<TH2D>(
			(std::string("h_BGOrec_shower_profile_cone_from_09_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("Shower Profile 0.9 < cosine < 1 - Reco Energy bin") + std::to_string(bidx + 1)).c_str(),
			14, 0, 13,
			10, 0, 1);

		h_BGOrec_layer_energy_ratio[bidx].resize(DAMPE_bgo_nLayers);
		h_BGOrec_layer_rms[bidx].resize(DAMPE_bgo_nLayers);
		h_BGOrec_energy_frac_1R[bidx].resize(DAMPE_bgo_nLayers);
		h_BGOrec_energy_frac_2R[bidx].resize(DAMPE_bgo_nLayers);
		h_BGOrec_energy_frac_3R[bidx].resize(DAMPE_bgo_nLayers);
		h_BGOrec_energy_frac_5R[bidx].resize(DAMPE_bgo_nLayers);

		for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
		{
			// h_BGOrec_layer_energy_ratio
			h_BGOrec_layer_energy_ratio[bidx][lIdx] = std::make_unique<TH1D>(
				(std::string("h_BGOrec_layer_energy_ratio_") + std::to_string(lIdx) + std::string("_recobin_") + std::to_string(bidx + 1)).c_str(),
				(std::string("Energy Ratio - BGO layer ") + std::to_string(lIdx) + std::string(" - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
				10, 0, 1);
			// h_BGOrec_layer_rms
			h_BGOrec_layer_rms[bidx][lIdx] = std::make_unique<TH1D>(
				(std::string("h_BGOrec_layer_rms_") + std::to_string(lIdx) + std::string("_recobin_") + std::to_string(bidx + 1)).c_str(),
				(std::string("RMS - BGO layer ") + std::to_string(lIdx) + std::string(" - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
				10, 0, 300);
			// h_BGOrec_energy_frac_1R
			h_BGOrec_energy_frac_1R[bidx][lIdx] = std::make_unique<TH1D>(
				(std::string("h_BGOrec_energy_frac_1R_") + std::to_string(lIdx) + std::string("_recobin_") + std::to_string(bidx + 1)).c_str(),
				(std::string("Energy fraction 1 Moliere Radius - BGO layer ") + std::to_string(lIdx) + std::string(" - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
				10, 0, 1);
			// h_BGOrec_energy_frac_2R
			h_BGOrec_energy_frac_2R[bidx][lIdx] = std::make_unique<TH1D>(
				(std::string("h_BGOrec_energy_frac_2R_") + std::to_string(lIdx) + std::string("_recobin_") + std::to_string(bidx + 1)).c_str(),
				(std::string("Energy fraction 2 Moliere Radius - BGO layer ") + std::to_string(lIdx) + std::string(" - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
				10, 0, 1);
			// h_BGOrec_energy_frac_3R
			h_BGOrec_energy_frac_3R[bidx][lIdx] = std::make_unique<TH1D>(
				(std::string("h_BGOrec_energy_frac_3R_") + std::to_string(lIdx) + std::string("_recobin_") + std::to_string(bidx + 1)).c_str(),
				(std::string("Energy fraction 3 Moliere Radius - BGO layer ") + std::to_string(lIdx) + std::string(" - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
				10, 0, 1);
			// h_BGOrec_energy_frac_5R
			h_BGOrec_energy_frac_5R[bidx][lIdx] = std::make_unique<TH1D>(
				(std::string("h_BGOrec_energy_frac_5R_") + std::to_string(lIdx) + std::string("_recobin_") + std::to_string(bidx + 1)).c_str(),
				(std::string("Energy fraction 5 Moliere Radius - BGO layer ") + std::to_string(lIdx) + std::string(" - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
				10, 0, 1);
		}

		h_BGOrec_sumrms[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_sumrms_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("sumRMS - Reco Energy bin ") + std::to_string(bidx + 1) + std::string("; sumRms[mm]")).c_str(),
			sumRms_bins.size() - 1, &(sumRms_bins[0]));
		h_BGOrec_sumrms_weighted[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_sumrms_weighted_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("weighted sumRMS - Reco Energy bin ") + std::to_string(bidx + 1) + std::string("; sumRms[mm]")).c_str(),
			sumRms_bins.size() - 1, &(sumRms_bins[0]));
		h_BGOrec_sumrms_cosine[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_sumrms_cosine_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("cosine sumRMS - Reco Energy bin ") + std::to_string(bidx + 1) + std::string("; sumRms[mm]")).c_str(),
			sumRms_bins.size() - 1, &(sumRms_bins[0]));
		h_BGOrec_fraclast[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_fraclast_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("Energy Fraction - BGO last layer - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
			100, 0, 1);
		h_BGOrec_fraclast_cosine[bidx] = std::make_unique<TH2D>(
			(std::string("h_BGOrec_fraclast_cosine_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("Energy Fraction - BGO last layer vs STK cosine Reco Energy bin ") + std::to_string(bidx + 1) + std::string("; cos(#theta); F_{last}")).c_str(),
			10, 0, 1,
			10, 0, 1);
		h_BGOrec_frac13[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_frac13_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("Energy Fraction - BGO 13th layer - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
			10, 0, 1);
		h_BGOrec_last_layer[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_last_layer_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("BGO last energy layer - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
			14, 0, 13);
		h_BGOrec_hits[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_hits_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("BGO hits - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
			10, 0, 1e+3);
		h_BGOrec_slopeX[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_slopeX_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("BGOrec Slope X - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
			10, -90, 90);
		h_BGOrec_slopeY[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_slopeY_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("BGOrec Slope Y - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
			10, -90, 90);
		h_BGOrec_interceptX[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_interceptX_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("BGOrec Intercept X - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
			50, -500, 500);
		h_BGOrec_interceptY[bidx] = std::make_unique<TH1D>(
			(std::string("h_BGOrec_interceptY_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("BGOrec Intercept Y - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
			50, -500, 500);
		h_BGOrec_topMap[bidx] = std::make_unique<TH2D>(
			(std::string("h_BGOrec_topMap_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("BGOreco TOP Map - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
			50, -500, 500,
			50, -500, 500);
		h_BGOrec_bottomMap[bidx] = std::make_unique<TH2D>(
			(std::string("h_BGOrec_bottomMap_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("BGOreco BOTTOM Map - Reco Energy bin ") + std::to_string(bidx + 1)).c_str(),
			50, -500, 500,
			50, -500, 500);
		sumRms_cosine[bidx] = std::make_unique<TH2D>(
			(std::string("sumRms_cosine_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("sumRms - cos(#theta) correlation - Reco Energy bin ") + std::to_string(bidx + 1) + std::string(" ; cos(#theta); sumRms [mm]")).c_str(),
			cosine_bins.size() - 1, &(cosine_bins[0]),
			sumRms_bins.size() - 1, &(sumRms_bins[0]));
		e_discrimination_last[bidx] = std::make_unique<TH2D>(
			(std::string("e_discrimination_last_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("Electron Discrimination - Reco Energy bin ") + std::to_string(bidx + 1) + std::string("; sumRms [mm]; F_{last}")).c_str(),
			sumRms_bins.size() - 1, &(sumRms_bins[0]),
			flast_binning.size() - 1, &(flast_binning[0]));
		e_discrimination[bidx] = std::make_unique<TH2D>(
			(std::string("e_discrimination_recobin_") + std::to_string(bidx + 1)).c_str(),
			(std::string("Electron Discrimination - Reco Energy bin ") + std::to_string(bidx + 1) + std::string("; sumRms [mm]; F")).c_str(),
			sumRms_bins.size() - 1, &(sumRms_bins[0]),
			flast_binning.size() - 1, &(flast_binning[0]));
	}

	sumRms_cosine_20_100 = std::make_unique<TH2D>(
		"sumRms_cosine_20_100",
		"sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_100_250 = std::make_unique<TH2D>(
		"sumRms_cosine_100_250", "sumRms - cos(#theta) correlation 100 GeV - 250 GeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_250_500 = std::make_unique<TH2D>(
		"sumRms_cosine_250_500",
		"sumRms - cos(#theta) correlation 250 GeV - 500 GeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_500_1000 = std::make_unique<TH2D>(
		"sumRms_cosine_500_1000",
		"sumRms - cos(#theta) correlation 500 GeV - 1 TeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_1000_3000 = std::make_unique<TH2D>(
		"sumRms_cosine_1000_3000",
		"sumRms - cos(#theta) correlation 1 TeV - 3 TeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_3000_10000 = std::make_unique<TH2D>(
		"sumRms_cosine_3000_10000",
		"sumRms - cos(#theta) correlation 3 TeV - 10 TeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_10000_20000 = std::make_unique<TH2D>(
		"sumRms_cosine_10000_20000",
		"sumRms - cos(#theta) correlation 10 TeV - 20 TeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));

	h_BGOrec_shower_profile_fit_energy = std::make_unique<TH1D>(
		"h_BGOrec_shower_profile_fit_energy",
		"Energy - BGO Shower Profile Fit",
		logEBins.size() - 1, &(logEBins[0]));
	h_BGOrec_shower_profile_fit_a = std::make_unique<TH1D>(
		"h_BGOrec_shower_profile_fit_a",
		"a parameter - BGO Shower Profile Fit",
		10, 0, 1);
	h_BGOrec_shower_profile_fit_b = std::make_unique<TH1D>(
		"h_BGOrec_shower_profile_fit_b",
		"b parameter - BGO Shower Profile Fit",
		10, 0, 1);
	h_BGOrec_shower_profile_energy_diff = std::make_unique<TH1D>(
		"h_BGOrec_shower_profile_energy_diff",
		"Fit Energy - BGO Correct energy - BGO Shower Profile Fit",
		10, -100, 100);
}

void histos::init_xtrl_histos()
{
	h_xtrl_energy_int = std::make_unique<TH1D>(
		"h_xtrl_energy_int",
		"Energy integrated XTRL distribution",
		xtrl_bins.size() - 1, &(xtrl_bins[0]));
	h_xtrl = std::make_unique<TH2D>(
		"h_xtrl",
		"XTRL energy Distribution; Corercted energy [GeV]; xtrl",
		logEBins.size() - 1, &(logEBins[0]),
		xtrl_bins.size() - 1, &(xtrl_bins[0]));
	for (unsigned int idx = 0; idx < logEBins.size() - 1; ++idx)
	{
		std::string bin_xtrl_name = "h_xtrl_bin_" + std::to_string(idx);
		h_xtrl_bin[idx] = std::make_unique<TH1D>(
			bin_xtrl_name.c_str(),
			"XTRL bin distribution; xtrl; counts",
			100, 0, 150);
	}
}

void histos::init_ep_histos()
{
	e_discrimination_last_20_100 = std::make_unique<TH2D>(
		"e_discrimination_last_20_100",
		"Electron Discrimination 20 GeV - 100 GeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_100_250 = std::make_unique<TH2D>(
		"e_discrimination_last_100_250",
		"Electron Discrimination 100 GeV - 250 GeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_250_500 = std::make_unique<TH2D>(
		"e_discrimination_last_250_500",
		"Electron Discrimination 250 GeV - 500 GeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_500_1000 = std::make_unique<TH2D>(
		"e_discrimination_last_500_1000",
		"Electron Discrimination 500 GeV - 1 TeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_1000_3000 = std::make_unique<TH2D>(
		"e_discrimination_last_1000_3000",
		"Electron Discrimination 1 TeV - 3 TeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_3000_10000 = std::make_unique<TH2D>(
		"e_discrimination_last_3000_10000",
		"Electron Discrimination 3 TeV - 10 TeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_10000_20000 = std::make_unique<TH2D>(
		"e_discrimination_last_10000_20000",
		"Electron Discrimination 10 TeV - 20 TeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));

	e_discrimination_20_100 = std::make_unique<TH2D>(
		"e_discrimination_20_100",
		"Electron Discrimination 20 GeV - 100 GeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_100_250 = std::make_unique<TH2D>(
		"e_discrimination_100_250",
		"Electron Discrimination 100 GeV - 250 GeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_250_500 = std::make_unique<TH2D>(
		"e_discrimination_250_500",
		"Electron Discrimination 250 GeV - 500 GeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_500_1000 = std::make_unique<TH2D>(
		"e_discrimination_500_1000",
		"Electron Discrimination 500 GeV - 1 TeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_1000_3000 = std::make_unique<TH2D>(
		"e_discrimination_1000_3000",
		"Electron Discrimination 1 TeV - 3 TeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_3000_10000 = std::make_unique<TH2D>(
		"e_discrimination_3000_10000",
		"Electron Discrimination 3 TeV - 10 TeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_10000_20000 = std::make_unique<TH2D>(
		"e_discrimination_10000_20000",
		"Electron Discrimination 10 TeV - 20 TeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
}

void histos::init_STK_histos()
{
	h_STK_costheta = std::make_unique<TH1D>(
		"h_STK_costheta",
		"STK costheta",
		10, 0, 1);
	h_STK_BGO_costheta_diff = std::make_unique<TH1D>(
		"h_STK_BGO_costheta_diff",
		"BGO reco - STK costheta",
		10, 0, 1);
}

void histos::init_psd_charge_histos()
{
	h_psd_chargeX = std::make_unique<TH1D>(
		"h_psd_chargeX",
		"Charge distribution X; X charge",
		10, 0, 100);
	h_psd_chargeY = std::make_unique<TH1D>(
		"h_psd_chargeY",
		"Charge distribution Y; Y charge",
		10, 0, 100);
	h_psd_charge2D = std::make_unique<TH2D>(
		"h_psd_charge2D",
		"PSD charge; X charge; Y charge",
		10, 0, 100,
		10, 0, 100);
	h_psd_charge = std::make_unique<TH1D>(
		"h_psd_charge",
		"Mean PSD charge; Mean charge",
		10, 0, 100);
	h_psd_selected_chargeX = std::make_unique<TH1D>(
		"h_psd_selected_chargeX",
		"Charge distribution X; X charge",
		10, 0, 100);
	h_psd_selected_chargeY = std::make_unique<TH1D>(
		"h_psd_selected_chargeY",
		"Charge distribution Y; Y charge",
		10, 0, 100);
	h_psd_selected_charge2D = std::make_unique<TH2D>(
		"h_psd_selected_charge2D",
		"PSD charge; X charge; Y charge",
		10, 0, 100,
		10, 0, 100);
	h_psd_selected_charge = std::make_unique<TH1D>(
		"h_psd_selected_charge",
		"Mean PSD charge; Mean charge",
		10, 0, 1000);
}

void histos::init_stk_charge_histos()
{
	h_stk_chargeX = std::make_unique<TH1D>(
		"h_stk_chargeX",
		"Charge distribution X; X charge",
		10, 0, 100);
	h_stk_chargeY = std::make_unique<TH1D>(
		"h_stk_chargeY",
		"Charge distribution Y; Y charge",
		10, 0, 100);
	h_stk_charge2D = std::make_unique<TH2D>(
		"h_stk_charge2D",
		"STK charge; X charge; Y charge",
		10, 0, 100,
		10, 0, 100);
	h_stk_charge = std::make_unique<TH1D>(
		"h_stk_charge",
		"Mean STK charge; Mean charge",
		10, 0, 100);
	h_stk_selected_chargeX = std::make_unique<TH1D>(
		"h_stk_selected_chargeX",
		"Charge distribution X; X charge",
		10, 0, 100);
	h_stk_selected_chargeY = std::make_unique<TH1D>(
		"h_stk_selected_chargeY",
		"Charge distribution Y; Y charge",
		10, 0, 100);
	h_stk_selected_charge2D = std::make_unique<TH2D>(
		"h_stk_selected_charge2D",
		"STK charge; X charge; Y charge",
		10, 0, 100,
		10, 0, 100);
	h_stk_selected_charge = std::make_unique<TH1D>(
		"h_stk_selected_charge",
		"Mean STK charge; Mean charge",
		10, 0, 100);
}

void histos::init_nud_histos()
{
	h_NUD_adc.resize(DAMPE_NUD_channels);

	for (int channel = 0; channel < DAMPE_NUD_channels; ++channel)
	{
		std::string h_name = "h_NUD_adc_" + std::to_string(channel);
		std::string h_title = "NUD ADC - channel " + std::to_string(channel);
		h_NUD_adc[channel] = std::make_unique<TH1D>(
			h_name.c_str(),
			h_title.c_str(),
			10, 0, 1e+3);
	}
	h_NUD_total_adc = std::make_unique<TH1D>(
		"h_NUD_total_adc",
		"NUD total ADC",
		100, 0, 1e+4);
	h_NUD_max_adc = std::make_unique<TH1D>(
		"h_NUD_max_adc",
		"NUD max ADC",
		10, 0, 1e+3);
	h_NUD_max_channel = std::make_unique<TH1D>(
		"h_NUD_max_channel",
		"NUD channel max ADC",
		4, 0, 3);
}

void histos::FillTrigger(
	const double energy,
	const double energy_w)
{
	h_trigger->Fill(energy, energy_w);
}

void histos::FillBGO(
	const double raw_energy,
	const double corr_energy,
	const double simu_energy_w,
	const int reco_bidx,
	const std::vector<double> *fracLayer,
	const std::vector<double> *eLayer,
	const std::vector<double> *rmsLayer,
	const std::vector<double> *energy_1R_radius,
	const std::vector<double> *energy_2R_radius,
	const std::vector<double> *energy_3R_radius,
	const std::vector<double> *energy_5R_radius,
	const double sumRms,
	const double fracLast,
	const double fracLast_13,
	const int lastBGOLayer,
	const int nBGOentries,
	const double BGOrec_slopeX,
	const double BGOrec_slopeY,
	const double BGOrec_interceptX,
	const double BGOrec_interceptY,
	const double BGOrec_costheta)
{
	const double _GeV = 0.001;

	h_BGOrec_energy->Fill(raw_energy * _GeV, simu_energy_w);
	h_BGOrec_corr_energy->Fill(corr_energy * _GeV, simu_energy_w);

	auto it_max = std::max_element(fracLayer->begin(), fracLayer->end());
	h_BGOrec_layer_max_energy_ratio[reco_bidx]->Fill(fracLayer->at(std::distance(fracLayer->begin(), it_max)), simu_energy_w);

	double weighted_sumRms = 0;
	double sum_energy_layer = 0;
	for (int lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
	{
		h_BGOrec_layer_energy_ratio[reco_bidx][lIdx]->Fill(fracLayer->at(lIdx), simu_energy_w);
		h_BGOrec_layer_rms[reco_bidx][lIdx]->Fill(rmsLayer->at(lIdx), simu_energy_w);
		h_BGOrec_energy_frac_1R[reco_bidx][lIdx]->Fill(energy_1R_radius->at(lIdx) / eLayer->at(lIdx), simu_energy_w);
		h_BGOrec_energy_frac_2R[reco_bidx][lIdx]->Fill(energy_2R_radius->at(lIdx) / eLayer->at(lIdx), simu_energy_w);
		h_BGOrec_energy_frac_3R[reco_bidx][lIdx]->Fill(energy_3R_radius->at(lIdx) / eLayer->at(lIdx), simu_energy_w);
		h_BGOrec_energy_frac_5R[reco_bidx][lIdx]->Fill(energy_5R_radius->at(lIdx) / eLayer->at(lIdx), simu_energy_w);
		h_BGOrec_shower_profile[reco_bidx]->Fill(lIdx, fracLayer->at(lIdx), simu_energy_w);
		if (BGOrec_costheta < .9)
			h_BGOrec_shower_profile_cosine_upto_09[reco_bidx]->Fill(lIdx, fracLayer->at(lIdx), simu_energy_w);
		else
			h_BGOrec_shower_profile_cone_from_09[reco_bidx]->Fill(lIdx, fracLayer->at(lIdx), simu_energy_w);

		weighted_sumRms += rmsLayer->at(lIdx) * eLayer->at(lIdx);
		sum_energy_layer += eLayer->at(lIdx);
	}

	weighted_sumRms /= raw_energy;

	h_BGOrec_sumrms[reco_bidx]->Fill(sumRms, simu_energy_w);
	h_BGOrec_sumrms_weighted[reco_bidx]->Fill(weighted_sumRms, simu_energy_w);
	h_BGOrec_sumrms_cosine[reco_bidx]->Fill(sumRms / BGOrec_costheta, simu_energy_w);
	h_BGOrec_layer_energy_diff->Fill(fabs(sum_energy_layer - raw_energy));
	h_BGOrec_fraclast[reco_bidx]->Fill(fracLast, simu_energy_w);
	h_BGOrec_fraclast_cosine[reco_bidx]->Fill(BGOrec_costheta, fracLast, simu_energy_w);
	h_BGOrec_frac13[reco_bidx]->Fill(fracLast_13, simu_energy_w);
	h_BGOrec_last_layer[reco_bidx]->Fill(lastBGOLayer, simu_energy_w);
	h_BGOrec_hits[reco_bidx]->Fill(nBGOentries, simu_energy_w);
	h_BGOrec_slopeX[reco_bidx]->Fill(BGOrec_slopeX, simu_energy_w);
	h_BGOrec_slopeY[reco_bidx]->Fill(BGOrec_slopeY, simu_energy_w);
	h_BGOrec_interceptX[reco_bidx]->Fill(BGOrec_interceptX, simu_energy_w);
	h_BGOrec_interceptY[reco_bidx]->Fill(BGOrec_interceptY, simu_energy_w);

	double BGOrec_topX = BGOrec_slopeX * BGO_TopZ + BGOrec_interceptX;
	double BGOrec_topY = BGOrec_slopeY * BGO_TopZ + BGOrec_interceptY;
	double BGOrec_bottomX = BGOrec_slopeX * BGO_BottomZ + BGOrec_interceptX;
	double BGOrec_bottomY = BGOrec_slopeY * BGO_BottomZ + BGOrec_interceptY;

	h_BGOrec_topMap[reco_bidx]->Fill(BGOrec_topX, BGOrec_topY, simu_energy_w);
	h_BGOrec_bottomMap[reco_bidx]->Fill(BGOrec_bottomX, BGOrec_bottomY, simu_energy_w);
	h_BGOrec_costheta->Fill(BGOrec_costheta, simu_energy_w);
}

void histos::FillBGOShowerFit(
	const std::vector<double> bgo_profile_fit_res,
	const double corr_energy,
	const double simu_energy_w)
{
	const double _GeV = 0.001;
	h_BGOrec_shower_profile_fit_energy->Fill(bgo_profile_fit_res[0] * _GeV, simu_energy_w);
	h_BGOrec_shower_profile_fit_b->Fill(bgo_profile_fit_res[1], simu_energy_w);
	h_BGOrec_shower_profile_fit_a->Fill(bgo_profile_fit_res[2], simu_energy_w);
	h_BGOrec_shower_profile_energy_diff->Fill((bgo_profile_fit_res[0] - corr_energy) * _GeV, simu_energy_w);
}

void histos::FillSumRmsCosine(
	const double energy,
	const double energy_w,
	const int reco_bidx,
	const double sum_rms,
	const double cosine)
{
	sumRms_cosine[reco_bidx]->Fill(cosine, sum_rms, energy_w);
	if (energy >= 20 && energy < 100)
		sumRms_cosine_20_100->Fill(cosine, sum_rms, energy_w);
	else if (energy >= 100 && energy < 250)
		sumRms_cosine_100_250->Fill(cosine, sum_rms, energy_w);
	else if (energy >= 250 && energy < 500)
		sumRms_cosine_250_500->Fill(cosine, sum_rms, energy_w);
	else if (energy >= 500 && energy < 1000)
		sumRms_cosine_500_1000->Fill(cosine, sum_rms, energy_w);
	else if (energy >= 1000 && energy < 3000)
		sumRms_cosine_1000_3000->Fill(cosine, sum_rms, energy_w);
	else if (energy >= 3000 && energy < 10000)
		sumRms_cosine_3000_10000->Fill(cosine, sum_rms, energy_w);
	else if (energy >= 10000 && energy < 20000)
		sumRms_cosine_10000_20000->Fill(cosine, sum_rms, energy_w);
}

void histos::FillSumRmsFLast(
	const double energy,
	const double energy_w,
	const int reco_bidx,
	const double sum_rms,
	const double last_energy_fraction,
	const double energy_fraction_13)
{
	if (last_energy_fraction != -1 && energy_fraction_13 != -1)
	{
		e_discrimination[reco_bidx]->Fill(sum_rms, last_energy_fraction, energy_w);
		e_discrimination_last[reco_bidx]->Fill(sum_rms, energy_fraction_13, energy_w);
		if (energy >= 20 && energy < 100)
		{
			e_discrimination_20_100->Fill(sum_rms, last_energy_fraction, energy_w);
			e_discrimination_last_20_100->Fill(sum_rms, energy_fraction_13, energy_w);
		}
		else if (energy >= 100 && energy < 250)
		{
			e_discrimination_100_250->Fill(sum_rms, last_energy_fraction, energy_w);
			e_discrimination_last_100_250->Fill(sum_rms, energy_fraction_13, energy_w);
		}
		else if (energy >= 250 && energy < 500)
		{
			e_discrimination_250_500->Fill(sum_rms, last_energy_fraction, energy_w);
			e_discrimination_last_250_500->Fill(sum_rms, energy_fraction_13, energy_w);
		}
		else if (energy >= 500 && energy < 1000)
		{
			e_discrimination_500_1000->Fill(sum_rms, last_energy_fraction, energy_w);
			e_discrimination_last_500_1000->Fill(sum_rms, energy_fraction_13, energy_w);
		}
		else if (energy >= 1000 && energy < 3000)
		{
			e_discrimination_1000_3000->Fill(sum_rms, last_energy_fraction, energy_w);
			e_discrimination_last_1000_3000->Fill(sum_rms, energy_fraction_13, energy_w);
		}
		else if (energy >= 3000 && energy < 10000)
		{
			e_discrimination_3000_10000->Fill(sum_rms, last_energy_fraction, energy_w);
			e_discrimination_last_3000_10000->Fill(sum_rms, energy_fraction_13, energy_w);
		}
		else if (energy >= 10000 && energy < 20000)
		{
			e_discrimination_10000_20000->Fill(sum_rms, last_energy_fraction, energy_w);
			e_discrimination_last_10000_20000->Fill(sum_rms, energy_fraction_13, energy_w);
		}
	}
}

void histos::FillGeometric(
	const double energy,
	const double energy_w)
{
	h_geometric_cut->Fill(energy, energy_w);
}

void histos::FillGeometricMaxLayer(
	const double energy,
	const double energy_w)
{
	h_geometric_maxElayer_cut->Fill(energy, energy_w);
}

void histos::FillGeometricMaxBar(
	const double energy,
	const double energy_w)
{
	h_geometric_maxBarLayer_cut->Fill(energy, energy_w);
}

void histos::FillGeometricBGOTrack(
	const double energy,
	const double energy_w)
{
	h_geometric_BGOTrackContainment_cut->Fill(energy, energy_w);
}

void histos::FillGeometricBGOFiducial(
	const double energy,
	const double energy_w)
{
	h_geometric_BGO_fiducial_cut->Fill(energy, energy_w);
}

void histos::FillGeometricAll(
	const double energy,
	const double energy_w)
{
	h_geometric_all_cut->Fill(energy, energy_w);
}

void histos::FillMaxLayer(
	const double energy,
	const double energy_w)
{
	h_maxElayer_cut->Fill(energy, energy_w);
}

void histos::FillMaxBar(
	const double energy,
	const double energy_w)
{
	h_maxBarLayer_cut->Fill(energy, energy_w);
}

void histos::FillBGOTrack(
	const double energy,
	const double energy_w)
{
	h_BGOTrackContainment_cut->Fill(energy, energy_w);
}

void histos::FillBGOFiducial(
	const double energy,
	const double energy_w)
{
	h_BGO_fiducial_cut->Fill(energy, energy_w);
}

void histos::FillBGOFiducialBarL13(
	const double energy,
	const double energy_w)
{
	h_BGOfiducial_nBarLayer13_cut->Fill(energy, energy_w);
}

void histos::FillBGOFiducialMaxRms(
	const double energy,
	const double energy_w)
{
	h_BGOfiducial_maxRms_cut->Fill(energy, energy_w);
}

void histos::FillBGOFiducialTrack(
	const double energy,
	const double energy_w)
{
	h_BGOfiducial_track_selection_cut->Fill(energy, energy_w);
}

void histos::FillBGOFiducialPsdStk(
	const double energy,
	const double energy_w)
{
	h_BGOfiducial_psd_stk_match_cut->Fill(energy, energy_w);
}

void histos::FillBGOFiducialPsdCharge(
	const double energy,
	const double energy_w)
{
	h_BGOfiducial_psd_charge_cut->Fill(energy, energy_w);
}

void histos::FillBGOFiducialStkCharge(
	const double energy,
	const double energy_w)
{
	h_BGOfiducial_stk_charge_cut->Fill(energy, energy_w);
}

void histos::FillBGOFiducialAll(
	const double energy,
	const double energy_w)
{
	h_BGOfiducial_all_cut->Fill(energy, energy_w);
}

void histos::FillBGOCosine(
	const double energy_w,
	const double BGOrec_cosine)
{
	h_BGOrec_costheta_presel->Fill(BGOrec_cosine, energy_w);
}

void histos::FillPsdCharge(
	const double energy_w,
	const double psd_charge_x,
	const double psd_charge_y,
	const double psd_charge,
	const bool selected)
{
	if (psd_charge_x != -999 && psd_charge_y != -999)
	{
		if (selected)
		{
			h_psd_selected_chargeX->Fill(psd_charge_x, energy_w);
			h_psd_selected_chargeY->Fill(psd_charge_y, energy_w);
			h_psd_selected_charge->Fill(psd_charge, energy_w);
			h_psd_selected_charge2D->Fill(psd_charge_x, psd_charge_y, energy_w);
		}
		else
		{
			h_psd_chargeX->Fill(psd_charge_x, energy_w);
			h_psd_chargeY->Fill(psd_charge_y, energy_w);
			h_psd_charge->Fill(psd_charge, energy_w);
			h_psd_charge2D->Fill(psd_charge_x, psd_charge_y, energy_w);
		}
	}
}

void histos::FillStkCharge(
	const double energy_w,
	const double stk_charge_x,
	const double stk_charge_y,
	const double stk_charge,
	const bool selected)
{
	if (stk_charge_x != -999 && stk_charge_y != -999)
	{
		if (selected)
		{
			h_stk_selected_chargeX->Fill(stk_charge_x, energy_w);
			h_stk_selected_chargeY->Fill(stk_charge_y, energy_w);
			h_stk_selected_charge->Fill(stk_charge, energy_w);
			h_stk_selected_charge2D->Fill(stk_charge_x, stk_charge_y, energy_w);
		}
		else
		{
			h_stk_chargeX->Fill(stk_charge_x, energy_w);
			h_stk_chargeY->Fill(stk_charge_y, energy_w);
			h_stk_charge->Fill(stk_charge, energy_w);
			h_stk_charge2D->Fill(stk_charge_x, stk_charge_y, energy_w);
		}
	}
}

void histos::FillStkCosine(
	const double energy_w,
	const double STK_cosine,
	const double BGOrec_cosine)
{
	h_STK_costheta->Fill(STK_cosine, energy_w);
	h_STK_BGO_costheta_diff->Fill(fabs(STK_cosine - BGOrec_cosine), energy_w);
}

void histos::FillClassifier(
	const double energy,
	const double energy_w,
	const double xtrl)
{
	h_xtrl_energy_int->Fill(xtrl, energy_w);
	h_xtrl->Fill(energy, xtrl, energy_w);
	auto xtrl_bin_idx = h_xtrl->GetXaxis()->FindBin(energy) - 1;
	if (xtrl_bin_idx >= 0 && (unsigned int)xtrl_bin_idx < logEBins.size() - 2)
		h_xtrl_bin[xtrl_bin_idx]->Fill(xtrl, energy_w);
}

void histos::FillNUD(
	const double energy_w,
	const std::vector<double> *nud_adc,
	const double nud_total_adc,
	const double nud_max_adc,
	const double nud_max_channel_id)
{
	for (int chIdx = 0; chIdx < DAMPE_NUD_channels; ++chIdx)
		h_NUD_adc[chIdx]->Fill(nud_adc->at(chIdx), energy_w);
	h_NUD_total_adc->Fill(nud_total_adc, energy_w);
	h_NUD_max_adc->Fill(nud_max_adc, energy_w);
	h_NUD_max_channel->Fill(nud_max_channel_id, energy_w);
}

void histos::FillAllCut(
	const double energy,
	const double energy_w)
{
	h_all_cut->Fill(energy, energy_w);
}

void histos::WriteCore(TFile *outfile)
{
	outfile->cd();

	auto trigger_dir = outfile->mkdir("Trigger");
	trigger_dir->cd();

	h_trigger->Write();

	auto preselection_dir = outfile->mkdir("Preselection");
	preselection_dir->cd();

	h_geometric_cut->Write();
	h_maxElayer_cut->Write();
	h_maxBarLayer_cut->Write();
	h_BGOTrackContainment_cut->Write();
	h_BGO_fiducial_cut->Write();
	h_all_cut->Write();

	h_geometric_maxElayer_cut->Write();
	h_geometric_maxBarLayer_cut->Write();
	h_geometric_BGOTrackContainment_cut->Write();
	h_geometric_BGO_fiducial_cut->Write();
	h_geometric_all_cut->Write();

	h_BGOfiducial_nBarLayer13_cut->Write();
	h_BGOfiducial_maxRms_cut->Write();
	h_BGOfiducial_track_selection_cut->Write();
	h_BGOfiducial_psd_stk_match_cut->Write();
	h_BGOfiducial_psd_charge_cut->Write();
	h_BGOfiducial_stk_charge_cut->Write();
	h_BGOfiducial_all_cut->Write();

	auto bgo_dir = outfile->mkdir("BGO");
	TDirectory *bin_tmp_dir = nullptr;
	std::string tmp_dir_name;
	for (unsigned int bidx = 0; bidx < (logEBins.size() - 1); ++bidx)
	{
		tmp_dir_name = std::string("BGO/recobin_") + std::to_string(bidx + 1);
		bin_tmp_dir = outfile->mkdir(tmp_dir_name.c_str());
		outfile->cd(tmp_dir_name.c_str());
		h_BGOrec_layer_max_energy_ratio[bidx]->Write();
		for (int lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
		{
			h_BGOrec_layer_energy_ratio[bidx][lIdx]->Write();
			h_BGOrec_layer_rms[bidx][lIdx]->Write();
			h_BGOrec_energy_frac_1R[bidx][lIdx]->Write();
			h_BGOrec_energy_frac_2R[bidx][lIdx]->Write();
			h_BGOrec_energy_frac_3R[bidx][lIdx]->Write();
			h_BGOrec_energy_frac_5R[bidx][lIdx]->Write();
		}

		h_BGOrec_sumrms[bidx]->Write();
		h_BGOrec_sumrms_weighted[bidx]->Write();
		h_BGOrec_sumrms_cosine[bidx]->Write();
		h_BGOrec_fraclast[bidx]->Write();
		h_BGOrec_fraclast_cosine[bidx]->Write();
		h_BGOrec_frac13[bidx]->Write();
		h_BGOrec_last_layer[bidx]->Write();
		h_BGOrec_hits[bidx]->Write();
		h_BGOrec_slopeX[bidx]->Write();
		h_BGOrec_slopeY[bidx]->Write();
		h_BGOrec_interceptX[bidx]->Write();
		h_BGOrec_interceptY[bidx]->Write();
		h_BGOrec_topMap[bidx]->Write();
		h_BGOrec_bottomMap[bidx]->Write();
		sumRms_cosine[bidx]->Write();
		e_discrimination_last[bidx]->Write();
		e_discrimination[bidx]->Write();
		h_BGOrec_shower_profile[bidx]->Write();
		h_BGOrec_shower_profile_cosine_upto_09[bidx]->Write();
		h_BGOrec_shower_profile_cone_from_09[bidx]->Write();
	}

	bgo_dir->cd();

	h_BGOrec_energy->Write();
	h_BGOrec_corr_energy->Write();
	h_BGOrec_layer_energy_diff->Write();
	h_BGOrec_costheta->Write();
	h_BGOrec_costheta_presel->Write();
	sumRms_cosine_20_100->Write();
	sumRms_cosine_100_250->Write();
	sumRms_cosine_250_500->Write();
	sumRms_cosine_500_1000->Write();
	sumRms_cosine_1000_3000->Write();
	sumRms_cosine_3000_10000->Write();
	sumRms_cosine_10000_20000->Write();

	e_discrimination_last_20_100->Write();
	e_discrimination_last_100_250->Write();
	e_discrimination_last_250_500->Write();
	e_discrimination_last_500_1000->Write();
	e_discrimination_last_1000_3000->Write();
	e_discrimination_last_3000_10000->Write();
	e_discrimination_last_10000_20000->Write();

	e_discrimination_20_100->Write();
	e_discrimination_100_250->Write();
	e_discrimination_250_500->Write();
	e_discrimination_500_1000->Write();
	e_discrimination_1000_3000->Write();
	e_discrimination_3000_10000->Write();
	e_discrimination_10000_20000->Write();

	h_BGOrec_shower_profile_fit_energy->Write();
	h_BGOrec_shower_profile_fit_a->Write();
	h_BGOrec_shower_profile_fit_b->Write();
	h_BGOrec_shower_profile_energy_diff->Write();

	auto psd_dir = outfile->mkdir("PSD");
	psd_dir->cd();

	h_psd_chargeX->Write();
	h_psd_chargeY->Write();
	h_psd_charge2D->Write();
	h_psd_charge->Write();
	h_psd_selected_chargeX->Write();
	h_psd_selected_chargeY->Write();
	h_psd_selected_charge2D->Write();
	h_psd_selected_charge->Write();

	auto stk_dir = outfile->mkdir("STK");
	stk_dir->cd();

	h_stk_chargeX->Write();
	h_stk_chargeY->Write();
	h_stk_charge2D->Write();
	h_stk_charge->Write();
	h_stk_selected_chargeX->Write();
	h_stk_selected_chargeY->Write();
	h_stk_selected_charge2D->Write();
	h_stk_selected_charge->Write();
	h_STK_costheta->Write();
	h_STK_BGO_costheta_diff->Write();

	auto nud_dir = outfile->mkdir("NUD");
	nud_dir->cd();

	for (auto &_elm : h_NUD_adc)
		_elm->Write();
	h_NUD_total_adc->Write();
	h_NUD_max_adc->Write();
	h_NUD_max_channel->Write();

	auto xtrl_dir = outfile->mkdir("xtrl");
	xtrl_dir->cd();

	h_xtrl_energy_int->Write();
	h_xtrl->Write();
	for (auto &_elm : h_xtrl_bin)
		_elm->Write();
}