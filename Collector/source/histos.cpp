#include "histos.h"
#include "binning.h"
#include "DAMPE_geo_structure.h"

histos::histos(std::vector<float> energy_bins)
{
	logEBins = energy_bins;
	init_preselection_histos();
	init_geometric_histos();
	init_BGOfiducial_histos();
	init_BGOlayer_histos();
	init_sumRms_cosine_histos();
	init_xtrl_histos();
	init_ep_histos();
	init_psd_charge_histos();
	init_stk_charge_histos();
	//init_time_histos();
}

void histos::init_preselection_histos()
{
	h_trigger = std::make_unique<TH1D>(
		"h_trigger",
		"Energy Distribution of the triggered particles; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_trigger->Sumw2();
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
	h_geometric_maxElayer_cut =std::make_unique<TH1D>(
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
	h_geo_BGOrec_slopeX = std::make_unique<TH1D>(
		"h_geo_BGOrec_slopeX",
		"BGOrec Slope X",
		1000, -90, 90);
	h_geo_BGOrec_slopeX->Sumw2();
	h_geo_BGOrec_slopeY = std::make_unique<TH1D>(
		"h_geo_BGOrec_slopeY",
		"BGOrec Slope Y",
		1000, -90, 90);
	h_geo_BGOrec_slopeY->Sumw2();
	h_geo_BGOrec_interceptX = std::make_unique<TH1D>(
		"h_geo_BGOrec_interceptX",
		"BGOrec Intercept X",
		500, -500, 500);
	h_geo_BGOrec_interceptX->Sumw2();
	h_geo_BGOrec_interceptY = std::make_unique<TH1D>(
		"h_geo_BGOrec_interceptY",
		"BGOrec Intercept Y",
		500, -500, 500);
	h_geo_BGOrec_interceptY->Sumw2();
	h_geo_BGOreco_topMap = std::make_unique<TH2D>(
		"h_geo_BGOreco_topMap",
		"BGOreco TOP Map",
		500, -500, 500,
		500, -500, 500);
	h_geo_BGOreco_topMap->Sumw2();
	h_geo_BGOreco_bottomMap = std::make_unique<TH2D>(
		"h_geo_BGOreco_bottomMap",
		"BGOreco BOTTOM Map",
		500, -500, 500,
		500, -500, 500);
	h_geo_BGOreco_bottomMap->Sumw2();
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
	h_BGOrec_energy = std::make_unique<TH1D>(
		"h_BGOrec_energy",
		"BGO Energy: Raw Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_BGOrec_energy->Sumw2();
	h_layer_max_energy_ratio = std::make_unique<TH1D>(
		"h_layer_max_energy_ratio",
		"Layer Energy Ratio",
		100, 0, 1);
	h_layer_max_energy_ratio->Sumw2();
}

void histos::init_BGOlayer_histos()
{
	h_layer_energy_ratio.resize(DAMPE_bgo_nLayers);
	std::string h_ratio_name;
	std::string h_ratio_title;
	for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
	{
		h_ratio_name = "h_layer_energy_ratio_" + std::to_string(lIdx);
		h_ratio_title = "Energy Ratio - BGO layer " + std::to_string(lIdx);
		h_layer_energy_ratio[lIdx] = std::make_unique<TH1D>(
			h_ratio_name.c_str(),
			h_ratio_title.c_str(),
			100, 0, 1);
		h_layer_energy_ratio[lIdx]->Sumw2();
	}
}

void histos::init_sumRms_cosine_histos()
{
	sumRms_cosine.resize(logEBins.size() - 1);
	cosine_bins = createLinearBinning(0, 1, 1e+2);
	sumRms_bins = createLogBinning(10, 2e+3, 1e+3);
	std::string histo_name;
	for (auto it = sumRms_cosine.begin(); it != sumRms_cosine.end(); ++it)
	{
		histo_name = "sumRms_cosine_" + to_string(std::distance(sumRms_cosine.begin(), it));
		(*it) = std::make_unique<TH2D>(
			histo_name.c_str(),
			"sumRms - cos(#theta) correlation; cos(#theta); sumRms [mm]",
			cosine_bins.size() - 1, &(cosine_bins[0]),
			sumRms_bins.size() - 1, &(sumRms_bins[0]));
		(*it)->Sumw2();
	}

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
}

void histos::init_xtrl_histos()
{
	xtrl_bins = createLinearBinning(0, 150, 1e+3);
	h_xtrl_bin.resize(logEBins.size() - 1);
	std::string bin_xtrl_name;
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
	for (unsigned int idx = 0; idx < logEBins.size(); ++idx)
	{
		bin_xtrl_name = "h_xtrl_bin_" + std::to_string(idx);
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

/*
void histos::init_time_histos()
{
	const int _tstamp_01012016 = 1451606400;
	const int _tstamp_31212020 = 1609459200;
	h_second = std::make_unique<TH1D>(
		"h_second", 
		"second", 
		100, time.start_second, time.end_second);
	h_second->Sumw2();
}
*/

void histos::evaluate_energy_ratio(
	const std::vector<double> fracLayer,
	const double energy_w)
{
	for (auto it = fracLayer.begin(); it != fracLayer.end(); ++it)
		h_layer_energy_ratio[std::distance(fracLayer.begin(), it)]->Fill(*it, energy_w);
	auto it_max = std::max_element(fracLayer.begin(), fracLayer.end());
	h_layer_max_energy_ratio->Fill(fracLayer[std::distance(fracLayer.begin(), it_max)], energy_w);
}

void histos::evaluate_top_bottom_position(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const double energy_w)
{
	double reco_topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
	double reco_topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

	double reco_bottomX = bgoRec_slope[0] * BGO_BottomZ + bgoRec_intercept[0];
	double reco_bottomY = bgoRec_slope[1] * BGO_BottomZ + bgoRec_intercept[1];

	// Fill slopes
	h_geo_BGOrec_slopeX->Fill(bgoRec_slope[0], energy_w);
	h_geo_BGOrec_slopeY->Fill(bgoRec_slope[1], energy_w);

	// Fill intercepts
	h_geo_BGOrec_interceptX->Fill(bgoRec_intercept[0], energy_w);
	h_geo_BGOrec_interceptY->Fill(bgoRec_intercept[1], energy_w);

	// Fill maps
	h_geo_BGOreco_topMap->Fill(reco_topX, reco_topY, energy_w);
	h_geo_BGOreco_bottomMap->Fill(reco_bottomX, reco_bottomY, energy_w);
}

void histos::fill_sumRms_cosine_histo(
	const double sumRMS,
	const double costheta,
	const double energy,
	const double energy_w)
{
	for (auto it = logEBins.begin(); it != logEBins.end() - 1; ++it)
		if (energy >= (*it) && energy < (*it + 1))
		{
			auto energy_idx = std::distance(logEBins.begin(), it);
			sumRms_cosine[energy_idx]->Fill(costheta, sumRMS, energy_w);
			break;
		}

	if (energy >= 20 && energy < 100)
		sumRms_cosine_20_100->Fill(costheta, sumRMS, energy_w);
	else if (energy >= 100 && energy < 250)
		sumRms_cosine_100_250->Fill(costheta, sumRMS, energy_w);
	else if (energy >= 250 && energy < 500)
		sumRms_cosine_250_500->Fill(costheta, sumRMS, energy_w);
	else if (energy >= 500 && energy < 1000)
		sumRms_cosine_500_1000->Fill(costheta, sumRMS, energy_w);
	else if (energy >= 1000 && energy < 3000)
		sumRms_cosine_1000_3000->Fill(costheta, sumRMS, energy_w);
	else
		sumRms_cosine_3000_10000->Fill(costheta, sumRMS, energy_w);
}

void histos::fill_XTRL_histo(
	const double xtrl,
	const double energy,
	const double energy_w)
{
	h_xtrl_energy_int->Fill(xtrl, energy_w);
	h_xtrl->Fill(energy, xtrl, energy_w);
	auto xtrl_bin_idx = h_xtrl->GetXaxis()->FindBin(energy) - 1;
	if (xtrl_bin_idx >= 0 && (unsigned int)xtrl_bin_idx < xtrl_bins.size())
		h_xtrl_bin[xtrl_bin_idx]->Fill(xtrl, energy_w);
}

void histos::fill_ep_histo(
	const double sumRMS,
	const double lastFracLayer,
	const double frac_layer_13,
	const double energy,
	const double energy_w)
{
	if (lastFracLayer != -1 && frac_layer_13 != -1)
	{
		e_discrimination->Fill(
			sumRMS,
			lastFracLayer,
			energy_w);
		e_discrimination_last->Fill(
			sumRMS,
			frac_layer_13,
			energy_w);

		if (energy >= 20 && energy < 100)
		{
			e_discrimination_20_100->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
			e_discrimination_last_20_100->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
		}
		else if (energy >= 100 && energy < 250)
		{
			e_discrimination_100_250->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
			e_discrimination_last_100_250->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
		}
		else if (energy >= 250 && energy < 500)
		{
			e_discrimination_250_500->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
			e_discrimination_last_250_500->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
		}
		else if (energy >= 500 && energy < 1000)
		{
			e_discrimination_500_1000->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
			e_discrimination_last_500_1000->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
		}
		else if (energy >= 1000 && energy < 3000)
		{
			e_discrimination_1000_3000->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
			e_discrimination_last_1000_3000->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
		}
		else
		{
			e_discrimination_3000_10000->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
			e_discrimination_last_3000_10000->Fill(
				sumRMS,
				lastFracLayer,
				energy_w);
		}
	}
}

void histos::Fill(
	const filter_output &output,
	const std::vector<double> &fracLayer,
	const std::vector<double> &bgoRec_slope,
	const std::vector<double> &bgoRec_intercept,
	const psd_charge &extracted_psd_charge,
	const stk_charge &extracted_stk_charge,
	const double sumRMS,
	const double costheta,
	const bgo_classifiers &classifier,
	const double lastFracLayer,
	const double frac_layer_13,
	const double energy,
	const double energy_w,
	const double simu_energy,
	const bool simu_evt)
{
	auto _GeV = 0.001;
	auto energy_gev = energy * _GeV;
	simu = simu_evt;
	double simu_energy_gev = -999;
	if (simu)
		simu_energy_gev = simu_energy * _GeV;
	double auto_energy_gev = -999;
	simu ? auto_energy_gev = simu_energy_gev : auto_energy_gev = energy_gev;

	// Check if the particle energy is outside SAA
	if (!output.evt_in_saa)
	{
		/*
		// Fill time histos for data events only
		if (!simu)
		{
			h_second->Fill(second);
			h_msecond->Fill(msecond);
		}
		*/
		// Check if the particle energy is within the energy range
		if (!output.out_energy_range)
		{
			// Check if the event has been triggered
			if (output.evt_triggered)
			{
				h_trigger->Fill(auto_energy_gev, energy_w);
				// Check if the event has been reconstructed succesfully
				if (output.correct_bgo_reco)
				{
					h_BGOrec_energy->Fill(energy_gev, energy_w);
					evaluate_energy_ratio(fracLayer, energy_w);

					if (output.geometric)
					{
						h_geometric_cut->Fill(auto_energy_gev, energy_w);
						evaluate_top_bottom_position(bgoRec_slope, bgoRec_intercept, energy_w);
						if (output.BGO_fiducial_maxElayer_cut)
							h_geometric_maxElayer_cut->Fill(
								auto_energy_gev,
								energy_w);
						if (output.BGO_fiducial_maxBarLayer_cut)
							h_geometric_maxBarLayer_cut->Fill(
								auto_energy_gev,
								energy_w);
						if (output.BGO_fiducial_BGOTrackContainment_cut)
							h_geometric_BGOTrackContainment_cut->Fill(
								auto_energy_gev,
								energy_w);
						if (output.BGO_fiducial)
							h_geometric_BGO_fiducial_cut->Fill(
								auto_energy_gev,
								energy_w);
						if (output.all_cut)
							h_geometric_all_cut->Fill(
								auto_energy_gev,
								energy_w);
					}

					if (output.BGO_fiducial_maxElayer_cut)
						h_maxElayer_cut->Fill(
							auto_energy_gev,
							energy_w);
					if (output.BGO_fiducial_maxBarLayer_cut)
						h_maxBarLayer_cut->Fill(
							auto_energy_gev,
							energy_w);
					if (output.BGO_fiducial_BGOTrackContainment_cut)
						h_BGOTrackContainment_cut->Fill(
							auto_energy_gev,
							energy_w);

					if (output.BGO_fiducial)
					{
						h_BGO_fiducial_cut->Fill(
							auto_energy_gev,
							energy_w);
						if (output.nBarLayer13_cut)
							h_BGOfiducial_nBarLayer13_cut->Fill(
								auto_energy_gev,
								energy_w);
						if (output.maxRms_cut)
							h_BGOfiducial_maxRms_cut->Fill(
								auto_energy_gev,
								energy_w);
						if (output.track_selection_cut)
							h_BGOfiducial_track_selection_cut->Fill(
								auto_energy_gev,
								energy_w);
						if (output.psd_stk_match_cut)
							h_BGOfiducial_psd_stk_match_cut->Fill(
								auto_energy_gev,
								energy_w);
						if (output.psd_charge_cut)
							h_BGOfiducial_psd_charge_cut->Fill(
								auto_energy_gev,
								energy_w);
						if (output.stk_charge_cut)
							h_BGOfiducial_stk_charge_cut->Fill(
								auto_energy_gev,
								energy_w);
						if (output.all_cut)
							h_BGOfiducial_all_cut->Fill(
								auto_energy_gev,
								energy_w);
					}

					if (output.psd_charge_measurement)
						if (extracted_psd_charge.chargeX != -999 && extracted_psd_charge.chargeY != -999)
						{
							auto mean_charge = 0.5 * (extracted_psd_charge.chargeX + extracted_psd_charge.chargeY);
							h_psd_chargeX->Fill(extracted_psd_charge.chargeX, energy_w);
							h_psd_chargeY->Fill(extracted_psd_charge.chargeY, energy_w);
							h_psd_charge->Fill(mean_charge, energy_w);
							h_psd_charge2D->Fill(extracted_psd_charge.chargeX, extracted_psd_charge.chargeY, energy_w);
						}
					if (output.psd_charge_cut)
						if (extracted_psd_charge.chargeX != -999 && extracted_psd_charge.chargeY != -999)
						{
							auto mean_charge = 0.5 * (extracted_psd_charge.chargeX + extracted_psd_charge.chargeY);
							h_psd_selected_chargeX->Fill(extracted_psd_charge.chargeX, energy_w);
							h_psd_selected_chargeY->Fill(extracted_psd_charge.chargeY, energy_w);
							h_psd_selected_charge->Fill(mean_charge, energy_w);
							h_psd_selected_charge2D->Fill(extracted_psd_charge.chargeX, extracted_psd_charge.chargeY, energy_w);
						}
					if (output.stk_charge_measurement)
						if (extracted_stk_charge.chargeX != -999 && extracted_stk_charge.chargeY != -999)
						{
							auto mean_charge = 0.5 * (extracted_stk_charge.chargeX + extracted_stk_charge.chargeY);
							h_stk_chargeX->Fill(extracted_stk_charge.chargeX, energy_w);
							h_stk_chargeY->Fill(extracted_stk_charge.chargeY, energy_w);
							h_stk_charge->Fill(mean_charge, energy_w);
							h_stk_charge2D->Fill(extracted_stk_charge.chargeX, extracted_stk_charge.chargeY, energy_w);
						}
					if (output.stk_charge_cut)
						if (extracted_stk_charge.chargeX != -999 && extracted_stk_charge.chargeY != -999)
						{
							auto mean_charge = 0.5 * (extracted_stk_charge.chargeX + extracted_stk_charge.chargeY);
							h_stk_selected_chargeX->Fill(extracted_stk_charge.chargeX, energy_w);
							h_stk_selected_chargeY->Fill(extracted_stk_charge.chargeY, energy_w);
							h_stk_selected_charge->Fill(mean_charge, energy_w);
							h_stk_selected_charge2D->Fill(extracted_stk_charge.chargeX, extracted_stk_charge.chargeY, energy_w);
						}

					if (output.all_cut)
					{
						h_all_cut->Fill(energy_gev, energy_w);
						fill_sumRms_cosine_histo(
							sumRMS,
							costheta,
							energy_gev,
							energy_w);
						fill_XTRL_histo(
							classifier.xtrl,
							energy_gev,
							energy_w);
						fill_ep_histo(
							sumRMS,
							lastFracLayer,
							frac_layer_13,
							energy_gev,
							energy_w);
					}
				}
			}
		}
	}
}

void histos::Write(TFile &outfile)
{
	outfile.cd();
	h_trigger->Write();
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

	auto bgorec_dir = outfile.mkdir("BGO_Rec");
	bgorec_dir->cd();

	h_geo_BGOrec_slopeX->Write();
	h_geo_BGOrec_slopeY->Write();
	h_geo_BGOrec_interceptX->Write();
	h_geo_BGOrec_interceptY->Write();
	h_geo_BGOreco_topMap->Write();
	h_geo_BGOreco_bottomMap->Write();

	auto bgoenergy_dir = outfile.mkdir("BGO_Energy");
	bgoenergy_dir->cd();

	h_BGOrec_energy->Write();
	h_layer_max_energy_ratio->Write();
	for (auto &_elm : h_layer_energy_ratio)
		_elm->Write();

	auto sumrms_dir = outfile.mkdir("SumRMS_cosine");
	sumrms_dir->cd();

	for (auto &_elm : sumRms_cosine)
		_elm->Write();
	sumRms_cosine_20_100->Write();
	sumRms_cosine_100_250->Write();
	sumRms_cosine_250_500->Write();
	sumRms_cosine_500_1000->Write();
	sumRms_cosine_1000_3000->Write();
	sumRms_cosine_3000_10000->Write();

	auto xtrl_dir = outfile.mkdir("xtrl");
	xtrl_dir->cd();

	h_xtrl_energy_int->Write();
	h_xtrl->Write();
	for (auto &_elm : h_xtrl_bin)
		_elm->Write();

	auto ep_dir = outfile.mkdir("ep_id");
	ep_dir->cd();

	e_discrimination_last->Write();
	e_discrimination_last_20_100->Write();
	e_discrimination_last_100_250->Write();
	e_discrimination_last_250_500->Write();
	e_discrimination_last_500_1000->Write();
	e_discrimination_last_1000_3000->Write();
	e_discrimination_last_3000_10000->Write();
	e_discrimination->Write();
	e_discrimination_20_100->Write();
	e_discrimination_100_250->Write();
	e_discrimination_250_500->Write();
	e_discrimination_500_1000->Write();
	e_discrimination_1000_3000->Write();
	e_discrimination_3000_10000->Write();

	auto psdcharge_dir = outfile.mkdir("PSDcharge");
	psdcharge_dir->cd();

	h_psd_chargeX->Write();
	h_psd_chargeY->Write();
	h_psd_charge->Write();
	h_psd_charge2D->Write();

	h_psd_selected_chargeX->Write();
	h_psd_selected_chargeY->Write();
	h_psd_selected_charge->Write();
	h_psd_selected_charge2D->Write();

	auto stkcharge_dir = outfile.mkdir("STKcharge");
	stkcharge_dir->cd();

	h_stk_chargeX->Write();
	h_stk_chargeY->Write();
	h_stk_charge->Write();
	h_stk_charge2D->Write();

	h_stk_selected_chargeX->Write();
	h_stk_selected_chargeY->Write();
	h_stk_selected_charge->Write();
	h_stk_selected_charge2D->Write();

	/*
	if (!simu)
	{
		auto time_dir = outfile.mkdir("Time");
		time_dir->cd();
		
		h_second->Write();
		h_msecond->Write();
	}
	*/
}