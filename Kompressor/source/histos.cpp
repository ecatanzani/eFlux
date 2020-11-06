#include "histos.h"
#include "binning.h"
#include "DAMPE_geo_structure.h"

#include "TDirectory.h"

histos::histos(std::vector<float> energy_bins)
{
	logEBins = energy_bins;
	init_trigger_histos();
	init_preselection_histos();
	init_geometric_histos();
	init_BGOfiducial_histos();
	init_BGO_histos();
	init_xtrl_histos();
	init_ep_histos();
	init_psd_charge_histos();
	init_stk_charge_histos();
	init_nud_histos();
}

void histos::init_trigger_histos()
{
	h_trigger = std::make_unique<TH1D>(
		"h_trigger",
		"Energy Distribution of the triggered particles; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_trigger->Sumw2();
}

void histos::init_preselection_histos()
{
	h_geometric_cut = std::make_unique<TH1D>(
		"h_geometric_cut",
		"Energy Distribution - geometric (trigger selection) cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_geometric_cut->Sumw2();
	h_maxElayer_cut = std::make_unique<TH1D>(
		"h_maxElayer_cut",
		"Energy Distribution - maxElayer cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_maxElayer_cut->Sumw2();
	h_maxBarLayer_cut = std::make_unique<TH1D>(
		"h_maxBarLayer_cut",
		"Energy Distribution - maxBarLayer cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_maxBarLayer_cut->Sumw2();
	h_BGOTrackContainment_cut = std::make_unique<TH1D>(
		"h_BGOTrackContainment_cut", "Energy Distribution - BGOTrackContainment cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOTrackContainment_cut->Sumw2();
	h_BGO_fiducial_cut = std::make_unique<TH1D>(
		"h_BGO_fiducial_cut", "Energy Distibution - BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGO_fiducial_cut->Sumw2();
	h_all_cut = std::make_unique<TH1D>(
		"h_all_cut", "Energy Distribution - All cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_all_cut->Sumw2();
}

void histos::init_geometric_histos()
{
	h_geometric_maxElayer_cut = std::make_unique<TH1D>(
		"h_geometric_maxElayer_cut",
		"Energy Distribution - maxElayer + geometric cut; Real Energy (GeV); counts",
		logEBins.size() - 1, &(logEBins[0]));
	h_geometric_maxElayer_cut->Sumw2();
	h_geometric_maxBarLayer_cut = std::make_unique<TH1D>(
		"h_geometric_maxBarLayer_cut",
		"Energy Distribution - maxBarLayer + geometric cut; Real Energy (GeV); counts",
		logEBins.size() - 1, &(logEBins[0]));
	h_geometric_maxBarLayer_cut->Sumw2();
	h_geometric_BGOTrackContainment_cut = std::make_unique<TH1D>(
		"h_geometric_BGOTrackContainment_cut",
		"Energy Distribution - BGOTrackContainment + geometric cut; Real Energy (GeV); counts",
		logEBins.size() - 1, &(logEBins[0]));
	h_geometric_BGOTrackContainment_cut->Sumw2();
	h_geometric_BGO_fiducial_cut = std::make_unique<TH1D>(
		"h_geometric_BGO_fiducial_cut",
		"Energy Distibution - BGO fiducial + geometric cut; Real Energy (GeV); counts",
		logEBins.size() - 1, &(logEBins[0]));
	h_geometric_BGO_fiducial_cut->Sumw2();
	h_geometric_all_cut = std::make_unique<TH1D>(
		"h_geometric_all_cut",
		"Energy Distribution - All + geometric cut; Real Energy (GeV); counts",
		logEBins.size() - 1, &(logEBins[0]));
	h_geometric_all_cut->Sumw2();
}

void histos::init_BGOfiducial_histos()
{
	h_BGOfiducial_nBarLayer13_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_nBarLayer13_cut",
		"Energy Distribution - nBarLayer13 + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_nBarLayer13_cut->Sumw2();
	h_BGOfiducial_maxRms_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_maxRms_cut",
		"Energy Distribution - maxRms  + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_maxRms_cut->Sumw2();
	h_BGOfiducial_track_selection_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_track_selection_cut",
		"Energy Distribution - track selection + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_track_selection_cut->Sumw2();
	h_BGOfiducial_psd_stk_match_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_psd_stk_match_cut",
		"Energy Distribution - psd-stk match + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_psd_stk_match_cut->Sumw2();
	h_BGOfiducial_psd_charge_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_psd_charge_cut",
		"Energy Distribution - psd charge + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_psd_charge_cut->Sumw2();
	h_BGOfiducial_stk_charge_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_stk_charge_cut",
		"Energy Distribution - stk charge + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_stk_charge_cut->Sumw2();
	h_BGOfiducial_all_cut = std::make_unique<TH1D>(
		"h_BGOfiducial_all_cut",
		"Energy Distribution - All + BGO fiducial cut; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOfiducial_all_cut->Sumw2();
}

void histos::init_BGO_histos()
{
	h_BGOrec_energy = std::make_unique<TH1D>(
		"h_BGOrec_energy",
		"BGO Energy: Raw Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOrec_energy->Sumw2();
	h_BGOrec_corr_energy = std::make_unique<TH1D>(
		"h_BGOrec_corr_energy",
		"BGO Corrected Energy: Raw Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOrec_corr_energy->Sumw2();
	h_BGOrec_layer_max_energy_ratio = std::make_unique<TH1D>(
		"h_BGOrec_layer_max_energy_ratio",
		"Max Layer Energy Ratio",
		100, 0, 1);
	h_BGOrec_layer_max_energy_ratio->Sumw2();

	h_BGOrec_layer_energy_ratio.resize(DAMPE_bgo_nLayers);
	h_BGOrec_layer_rms.resize(DAMPE_bgo_nLayers);
	h_BGOrec_energy_frac_1R.resize(DAMPE_bgo_nLayers);
	h_BGOrec_energy_frac_2R.resize(DAMPE_bgo_nLayers);
	h_BGOrec_energy_frac_3R.resize(DAMPE_bgo_nLayers);
	h_BGOrec_energy_frac_5R.resize(DAMPE_bgo_nLayers);

	for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
	{
		// h_BGOrec_layer_energy_ratio
		std::string h_name = "h_BGOrec_layer_energy_ratio_" + std::to_string(lIdx);
		std::string h_title = "Energy Ratio - BGO layer " + std::to_string(lIdx);
		h_BGOrec_layer_energy_ratio[lIdx] = std::make_unique<TH1D>(
			h_name.c_str(),
			h_title.c_str(),
			100, 0, 1);
		h_BGOrec_layer_energy_ratio[lIdx]->Sumw2();
		// h_BGOrec_layer_rms
		h_name = "h_BGOrec_layer_rms_" + std::to_string(lIdx);
		h_title = "RMS - BGO layer " + std::to_string(lIdx);
		h_BGOrec_layer_rms[lIdx] = std::make_unique<TH1D>(
			h_name.c_str(),
			h_title.c_str(),
			1000, 0, 300);
		h_BGOrec_layer_rms[lIdx]->Sumw2();
		// h_BGOrec_energy_frac_1R
		h_name = "h_BGOrec_energy_frac_1R_" + std::to_string(lIdx);
		h_title = "Energy fraction 1 Moliere Radius - BGO layer " + std::to_string(lIdx);
		h_BGOrec_energy_frac_1R[lIdx] = std::make_unique<TH1D>(
			h_name.c_str(),
			h_title.c_str(),
			100, 0, 1);
		h_BGOrec_energy_frac_1R[lIdx]->Sumw2();
		// h_BGOrec_energy_frac_2R
		h_name = "h_BGOrec_energy_frac_2R_" + std::to_string(lIdx);
		h_title = "Energy fraction 2 Moliere Radius - BGO layer " + std::to_string(lIdx);
		h_BGOrec_energy_frac_2R[lIdx] = std::make_unique<TH1D>(
			h_name.c_str(),
			h_title.c_str(),
			100, 0, 1);
		h_BGOrec_energy_frac_2R[lIdx]->Sumw2();
		// h_BGOrec_energy_frac_3R
		h_name = "h_BGOrec_energy_frac_3R_" + std::to_string(lIdx);
		h_title = "Energy fraction 3 Moliere Radius - BGO layer " + std::to_string(lIdx);
		h_BGOrec_energy_frac_3R[lIdx] = std::make_unique<TH1D>(
			h_name.c_str(),
			h_title.c_str(),
			100, 0, 1);
		h_BGOrec_energy_frac_3R[lIdx]->Sumw2();
		// h_BGOrec_energy_frac_5R
		h_name = "h_BGOrec_energy_frac_5R_" + std::to_string(lIdx);
		h_title = "Energy fraction 5 Moliere Radius - BGO layer " + std::to_string(lIdx);
		h_BGOrec_energy_frac_5R[lIdx] = std::make_unique<TH1D>(
			h_name.c_str(),
			h_title.c_str(),
			100, 0, 1);
		h_BGOrec_energy_frac_5R[lIdx]->Sumw2();
	}

	sumRms_bins = createLogBinning(10, 2e+3, 1e+3);

	h_BGOrec_sumrms = std::make_unique<TH1D>(
		"h_BGOrec_sumrms",
		"sumRMS; sumRms[mm]",
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	h_BGOrec_sumrms->Sumw2();
	h_BGOrec_fraclast = std::make_unique<TH1D>(
		"h_BGOrec_fraclast",
		"Energy Fraction - BGO last layer",
		100, 0, 1);
	h_BGOrec_fraclast->Sumw2();
	h_BGOrec_frac13 = std::make_unique<TH1D>(
		"h_BGOrec_frac13",
		"Energy Fraction - BGO 13th layer",
		100, 0, 1);
	h_BGOrec_frac13->Sumw2();
	h_BGOrec_last_layer = std::make_unique<TH1D>(
		"h_BGOrec_last_layer",
		"BGO last energy layer",
		14, 0, 13);
	h_BGOrec_last_layer->Sumw2();
	h_BGOrec_hits = std::make_unique<TH1D>(
		"h_BGOrec_hits",
		"BGO hits",
		1e+2, 0, 1e+3);
	h_BGOrec_hits->Sumw2();
	h_BGOrec_slopeX = std::make_unique<TH1D>(
		"h_BGOrec_slopeX",
		"BGOrec Slope X",
		1000, -90, 90);
	h_BGOrec_slopeX->Sumw2();
	h_BGOrec_slopeY = std::make_unique<TH1D>(
		"h_BGOrec_slopeY",
		"BGOrec Slope Y",
		1000, -90, 90);
	h_BGOrec_slopeY->Sumw2();
	h_BGOrec_interceptX = std::make_unique<TH1D>(
		"h_BGOrec_interceptX",
		"BGOrec Intercept X",
		500, -500, 500);
	h_BGOrec_interceptX->Sumw2();
	h_BGOrec_interceptY = std::make_unique<TH1D>(
		"h_BGOrec_interceptY",
		"BGOrec Intercept Y",
		500, -500, 500);
	h_BGOrec_interceptY->Sumw2();
	h_BGOrec_topMap = std::make_unique<TH2D>(
		"h_BGOrec_topMap",
		"BGOreco TOP Map",
		500, -500, 500,
		500, -500, 500);
	h_BGOrec_topMap->Sumw2();
	h_BGOrec_bottomMap = std::make_unique<TH2D>(
		"h_BGOrec_bottomMap",
		"BGOreco BOTTOM Map",
		500, -500, 500,
		500, -500, 500);
	h_BGOrec_bottomMap->Sumw2();

	cosine_bins = createLinearBinning(0, 1, 1e+2);

	sumRms_cosine_20_100 = std::make_unique<TH2D>(
		"sumRms_cosine_20_100",
		"sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_20_100->Sumw2();
	sumRms_cosine_100_250 = std::make_unique<TH2D>(
		"sumRms_cosine_100_250", "sumRms - cos(#theta) correlation 100 GeV - 250 GeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_100_250->Sumw2();
	sumRms_cosine_250_500 = std::make_unique<TH2D>(
		"sumRms_cosine_250_500",
		"sumRms - cos(#theta) correlation 250 GeV - 500 GeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_250_500->Sumw2();
	sumRms_cosine_500_1000 = std::make_unique<TH2D>(
		"sumRms_cosine_500_1000",
		"sumRms - cos(#theta) correlation 500 GeV - 1 TeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_500_1000->Sumw2();
	sumRms_cosine_1000_3000 = std::make_unique<TH2D>(
		"sumRms_cosine_1000_3000",
		"sumRms - cos(#theta) correlation 1 TeV - 3 TeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_1000_3000->Sumw2();
	sumRms_cosine_3000_10000 = std::make_unique<TH2D>(
		"sumRms_cosine_3000_10000",
		"sumRms - cos(#theta) correlation 3 TeV - 10 TeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_3000_10000->Sumw2();
	sumRms_cosine_10000_20000 = std::make_unique<TH2D>(
		"sumRms_cosine_10000_20000",
		"sumRms - cos(#theta) correlation 10 TeV - 20 TeV; cos(#theta); sumRms [mm]",
		cosine_bins.size() - 1, &(cosine_bins[0]),
		sumRms_bins.size() - 1, &(sumRms_bins[0]));
	sumRms_cosine_10000_20000->Sumw2();
}

void histos::init_xtrl_histos()
{
	xtrl_bins = createLinearBinning(0, 150, 1e+3);
	h_xtrl_bin.resize(logEBins.size() - 1);
	h_xtrl_energy_int = std::make_unique<TH1D>(
		"h_xtrl_energy_int",
		"Energy integrated XTRL distribution",
		xtrl_bins.size() - 1, &(xtrl_bins[0]));
	h_xtrl_energy_int->Sumw2();
	h_xtrl = std::make_unique<TH2D>(
		"h_xtrl",
		"XTRL energy Distribution; Corercted energy [GeV]; xtrl",
		logEBins.size() - 1, &(logEBins[0]),
		xtrl_bins.size() - 1, &(xtrl_bins[0]));
	h_xtrl->Sumw2();
	for (unsigned int idx = 0; idx < logEBins.size() - 1; ++idx)
	{
		std::string bin_xtrl_name = "h_xtrl_bin_" + std::to_string(idx);
		h_xtrl_bin[idx] = std::make_unique<TH1D>(
			bin_xtrl_name.c_str(),
			"XTRL bin distribution; xtrl; counts",
			100, 0, 150);
		h_xtrl_bin[idx]->Sumw2();
	}
}

void histos::init_ep_histos()
{
	flast_binning = createLogBinning(1e-5, 2e-1, 1e+3);
	e_discrimination_last = std::make_unique<TH2D>(
		"e_discrimination_last",
		"Electron Discrimination; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last->Sumw2();
	e_discrimination_last_20_100 = std::make_unique<TH2D>(
		"e_discrimination_last_20_100",
		"Electron Discrimination 20 GeV - 100 GeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_20_100->Sumw2();
	e_discrimination_last_100_250 = std::make_unique<TH2D>(
		"e_discrimination_last_100_250",
		"Electron Discrimination 100 GeV - 250 GeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_100_250->Sumw2();
	e_discrimination_last_250_500 = std::make_unique<TH2D>(
		"e_discrimination_last_250_500",
		"Electron Discrimination 250 GeV - 500 GeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_250_500->Sumw2();
	e_discrimination_last_500_1000 = std::make_unique<TH2D>(
		"e_discrimination_last_500_1000",
		"Electron Discrimination 500 GeV - 1 TeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_500_1000->Sumw2();
	e_discrimination_last_1000_3000 = std::make_unique<TH2D>(
		"e_discrimination_last_1000_3000",
		"Electron Discrimination 1 TeV - 3 TeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_1000_3000->Sumw2();
	e_discrimination_last_3000_10000 = std::make_unique<TH2D>(
		"e_discrimination_last_3000_10000",
		"Electron Discrimination 3 TeV - 10 TeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_3000_10000->Sumw2();
	e_discrimination_last_10000_20000 = std::make_unique<TH2D>(
		"e_discrimination_last_10000_20000",
		"Electron Discrimination 10 TeV - 20 TeV; sumRms [mm]; F_{last}",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_last_10000_20000->Sumw2();

	e_discrimination = std::make_unique<TH2D>(
		"e_discrimination",
		"Electron Discrimination; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination->Sumw2();
	e_discrimination_20_100 = std::make_unique<TH2D>(
		"e_discrimination_20_100",
		"Electron Discrimination 20 GeV - 100 GeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_20_100->Sumw2();
	e_discrimination_100_250 = std::make_unique<TH2D>(
		"e_discrimination_100_250",
		"Electron Discrimination 100 GeV - 250 GeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_100_250->Sumw2();
	e_discrimination_250_500 = std::make_unique<TH2D>(
		"e_discrimination_250_500",
		"Electron Discrimination 250 GeV - 500 GeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_250_500->Sumw2();
	e_discrimination_500_1000 = std::make_unique<TH2D>(
		"e_discrimination_500_1000",
		"Electron Discrimination 500 GeV - 1 TeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_500_1000->Sumw2();
	e_discrimination_1000_3000 = std::make_unique<TH2D>(
		"e_discrimination_1000_3000",
		"Electron Discrimination 1 TeV - 3 TeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_1000_3000->Sumw2();
	e_discrimination_3000_10000 = std::make_unique<TH2D>(
		"e_discrimination_3000_10000",
		"Electron Discrimination 3 TeV - 10 TeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_3000_10000->Sumw2();
	e_discrimination_10000_20000 = std::make_unique<TH2D>(
		"e_discrimination_10000_20000",
		"Electron Discrimination 10 TeV - 20 TeV; sumRms [mm]; F",
		sumRms_bins.size() - 1, &(sumRms_bins[0]),
		flast_binning.size() - 1, &(flast_binning[0]));
	e_discrimination_10000_20000->Sumw2();
}

void histos::init_psd_charge_histos()
{
	h_psd_chargeX = std::make_unique<TH1D>(
		"h_psd_chargeX",
		"Charge distribution X; X charge",
		1000, 0, 100);
	h_psd_chargeX->Sumw2();
	h_psd_chargeY = std::make_unique<TH1D>(
		"h_psd_chargeY",
		"Charge distribution Y; Y charge",
		1000, 0, 100);
	h_psd_chargeY->Sumw2();
	h_psd_charge2D = std::make_unique<TH2D>(
		"h_psd_charge2D",
		"PSD charge; X charge; Y charge",
		1000, 0, 100,
		1000, 0, 100);
	h_psd_charge2D->Sumw2();
	h_psd_charge = std::make_unique<TH1D>(
		"h_psd_charge",
		"Mean PSD charge; Mean charge",
		1000, 0, 100);
	h_psd_charge->Sumw2();
	h_psd_selected_chargeX = std::make_unique<TH1D>(
		"h_psd_selected_chargeX",
		"Charge distribution X; X charge",
		1000, 0, 100);
	h_psd_selected_chargeX->Sumw2();
	h_psd_selected_chargeY = std::make_unique<TH1D>(
		"h_psd_selected_chargeY",
		"Charge distribution Y; Y charge",
		1000, 0, 100);
	h_psd_selected_chargeY->Sumw2();
	h_psd_selected_charge2D = std::make_unique<TH2D>(
		"h_psd_selected_charge2D",
		"PSD charge; X charge; Y charge",
		1000, 0, 100,
		1000, 0, 100);
	h_psd_selected_charge2D->Sumw2();
	h_psd_selected_charge = std::make_unique<TH1D>(
		"h_psd_selected_charge",
		"Mean PSD charge; Mean charge",
		1000, 0, 1000);
	h_psd_selected_charge->Sumw2();
}

void histos::init_stk_charge_histos()
{
	h_stk_chargeX = std::make_unique<TH1D>(
		"h_stk_chargeX",
		"Charge distribution X; X charge",
		1000, 0, 100);
	h_stk_chargeX->Sumw2();
	h_stk_chargeY = std::make_unique<TH1D>(
		"h_stk_chargeY",
		"Charge distribution Y; Y charge",
		1000, 0, 100);
	h_stk_chargeY->Sumw2();
	h_stk_charge2D = std::make_unique<TH2D>(
		"h_stk_charge2D",
		"STK charge; X charge; Y charge",
		1000, 0, 100,
		1000, 0, 100);
	h_stk_charge2D->Sumw2();
	h_stk_charge = std::make_unique<TH1D>(
		"h_stk_charge",
		"Mean STK charge; Mean charge",
		1000, 0, 100);
	h_stk_charge->Sumw2();
	h_stk_selected_chargeX = std::make_unique<TH1D>(
		"h_stk_selected_chargeX",
		"Charge distribution X; X charge",
		1000, 0, 100);
	h_stk_selected_chargeX->Sumw2();
	h_stk_selected_chargeY = std::make_unique<TH1D>(
		"h_stk_selected_chargeY",
		"Charge distribution Y; Y charge",
		1000, 0, 100);
	h_stk_selected_chargeY->Sumw2();
	h_stk_selected_charge2D = std::make_unique<TH2D>(
		"h_stk_selected_charge2D",
		"STK charge; X charge; Y charge",
		1000, 0, 100,
		1000, 0, 100);
	h_stk_selected_charge2D->Sumw2();
	h_stk_selected_charge = std::make_unique<TH1D>(
		"h_stk_selected_charge",
		"Mean STK charge; Mean charge",
		1000, 0, 100);
	h_stk_selected_charge->Sumw2();
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
			100, 0, 1e+3);
		h_NUD_adc[channel]->Sumw2();
	}
	h_NUD_total_adc = std::make_unique<TH1D>(
		"h_NUD_total_adc",
		"NUD total ADC",
		1e+3, 0, 1e+4);
	h_NUD_total_adc->Sumw2();
	h_NUD_max_adc = std::make_unique<TH1D>(
		"h_NUD_max_adc",
		"NUD max ADC",
		10, 0, 1e+3);
	h_NUD_max_adc->Sumw2();
	h_NUD_max_channel = std::make_unique<TH1D>(
		"h_NUD_max_channel",
		"NUD channel max ADC",
		4, 0, 3);
	h_NUD_max_channel->Sumw2();
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
	const double BGOrec_interceptY)
{
	h_BGOrec_energy->Fill(raw_energy, simu_energy_w);
	h_BGOrec_corr_energy->Fill(corr_energy, simu_energy_w);

	auto it_max = std::max_element(fracLayer->begin(), fracLayer->end());
	h_BGOrec_layer_max_energy_ratio->Fill(fracLayer->at(std::distance(fracLayer->begin(), it_max)), simu_energy_w);

	for (int lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
	{
		h_BGOrec_layer_energy_ratio[lIdx]->Fill(fracLayer->at(lIdx), simu_energy_w);
		h_BGOrec_layer_rms[lIdx]->Fill(rmsLayer->at(lIdx), simu_energy_w);
		h_BGOrec_energy_frac_1R[lIdx]->Fill(energy_1R_radius->at(lIdx) / eLayer->at(lIdx), simu_energy_w);
		h_BGOrec_energy_frac_2R[lIdx]->Fill(energy_2R_radius->at(lIdx) / eLayer->at(lIdx), simu_energy_w);
		h_BGOrec_energy_frac_3R[lIdx]->Fill(energy_3R_radius->at(lIdx) / eLayer->at(lIdx), simu_energy_w);
		h_BGOrec_energy_frac_5R[lIdx]->Fill(energy_5R_radius->at(lIdx) / eLayer->at(lIdx), simu_energy_w);
	}

	h_BGOrec_sumrms->Fill(sumRms, simu_energy_w);
	h_BGOrec_fraclast->Fill(fracLast, simu_energy_w);
	h_BGOrec_frac13->Fill(fracLast_13, simu_energy_w);
	h_BGOrec_last_layer->Fill(lastBGOLayer, simu_energy_w);
	h_BGOrec_hits->Fill(nBGOentries, simu_energy_w);
	h_BGOrec_slopeX->Fill(BGOrec_slopeX, simu_energy_w);
	h_BGOrec_slopeY->Fill(BGOrec_slopeY, simu_energy_w);
	h_BGOrec_interceptX->Fill(BGOrec_interceptX, simu_energy_w);
	h_BGOrec_interceptY->Fill(BGOrec_interceptY, simu_energy_w);

	double BGOrec_topX = BGOrec_slopeX * BGO_TopZ + BGOrec_interceptX;
	double BGOrec_topY = BGOrec_slopeY * BGO_TopZ + BGOrec_interceptY;
	double BGOrec_bottomX = BGOrec_slopeX * BGO_BottomZ + BGOrec_interceptX;
	double BGOrec_bottomY = BGOrec_slopeY * BGO_BottomZ + BGOrec_interceptY;

	h_BGOrec_topMap->Fill(BGOrec_topX, BGOrec_topY, simu_energy_w);
	h_BGOrec_bottomMap->Fill(BGOrec_bottomX, BGOrec_bottomY, simu_energy_w);
}

void histos::FillSumRmsCosine(
	const double energy,
	const double energy_w,
	const double sum_rms,
	const double cosine)
{
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
	const double sum_rms,
	const double last_energy_fraction,
	const double energy_fraction_13)
{
	if (last_energy_fraction != -1 && energy_fraction_13 != -1)
	{
		e_discrimination->Fill(sum_rms, last_energy_fraction, energy_w);
		e_discrimination_last->Fill(sum_rms, energy_fraction_13, energy_w);
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
	bgo_dir->cd();

	h_BGOrec_energy->Write();
	h_BGOrec_corr_energy->Write();
	h_BGOrec_layer_max_energy_ratio->Write();

	for (int lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
	{
		h_BGOrec_layer_energy_ratio[lIdx]->Write();
		h_BGOrec_layer_rms[lIdx]->Write();
		h_BGOrec_energy_frac_1R[lIdx]->Write();
		h_BGOrec_energy_frac_2R[lIdx]->Write();
		h_BGOrec_energy_frac_3R[lIdx]->Write();
		h_BGOrec_energy_frac_5R[lIdx]->Write();
	}

	h_BGOrec_sumrms->Write();
	h_BGOrec_fraclast->Write();
	h_BGOrec_frac13->Write();
	h_BGOrec_last_layer->Write();
	h_BGOrec_hits->Write();
	h_BGOrec_slopeX->Write();
	h_BGOrec_slopeY->Write();
	h_BGOrec_interceptX->Write();
	h_BGOrec_interceptY->Write();
	h_BGOrec_topMap->Write();
	h_BGOrec_bottomMap->Write();

	sumRms_cosine_20_100->Write();
	sumRms_cosine_100_250->Write();
	sumRms_cosine_250_500->Write();
	sumRms_cosine_500_1000->Write();
	sumRms_cosine_1000_3000->Write();
	sumRms_cosine_3000_10000->Write();
	sumRms_cosine_10000_20000->Write();

	e_discrimination_last->Write();
	e_discrimination_last_20_100->Write();
	e_discrimination_last_100_250->Write();
	e_discrimination_last_250_500->Write();
	e_discrimination_last_500_1000->Write();
	e_discrimination_last_1000_3000->Write();
	e_discrimination_last_3000_10000->Write();
	e_discrimination_last_10000_20000->Write();

	e_discrimination->Write();
	e_discrimination_20_100->Write();
	e_discrimination_100_250->Write();
	e_discrimination_250_500->Write();
	e_discrimination_500_1000->Write();
	e_discrimination_1000_3000->Write();
	e_discrimination_3000_10000->Write();
	e_discrimination_10000_20000->Write();

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

	auto nud_dir = outfile->mkdir("NUD");
	nud_dir->cd();

	for (auto& _elm : h_NUD_adc)
		_elm->Write();
	h_NUD_total_adc->Write();
	h_NUD_max_adc->Write();
	h_NUD_max_channel->Write();

	auto xtrl_dir = outfile->mkdir("xtrl");
	xtrl_dir->cd();

	h_xtrl_energy_int->Write();
	h_xtrl->Write();
	for (auto& _elm : h_xtrl_bin)
		_elm->Write();
}