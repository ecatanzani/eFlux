#include "wtsydp.h"
#include "flux.h"
#include "data_cuts.h"

#include "TH2D.h"
#include "TGraphErrors.h"
#include "TEfficiency.h"

void generateDataFinalGraph(
	const bool verbose,
	const bool pedantic,
	const std::string outputPath,
	const std::string complete_histo_path,
	const std::string wd)
{
	if (verbose)
		std::cout << "\n*** Extractor DATA facility ***\n\n";
	
	// Create flux cuts struct
	cuts_conf flux_cuts;
	// Create active cuts struct
	data_active_cuts active_cuts;

	// Load structs reading config file
	auto logEBins = load_flux_struct(
		flux_cuts, 
		active_cuts, 
		wd);
	
	// Open input complete histos
	TFile inHisto(complete_histo_path.c_str(), "READ");
	if (!inHisto.IsOpen())
	{
		std::cerr << "Error reading input complete histo: " << complete_histo_path;
		exit(123);
	}

	auto h_trigger = static_cast<TH1D *>(inHisto.Get("h_trigger"));
	auto h_geometric_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_cut"));
	auto h_maxElayer_cut = static_cast<TH1D *>(inHisto.Get("h_maxElayer_cut"));
	auto h_maxBarLayer_cut = static_cast<TH1D *>(inHisto.Get("h_maxBarLayer_cut"));
	auto h_BGOTrackContainment_cut = static_cast<TH1D *>(inHisto.Get("h_BGOTrackContainment_cut"));
	auto h_BGO_fiducial_cut = static_cast<TH1D *>(inHisto.Get("h_BGO_fiducial_cut"));
	auto h_all_cut = static_cast<TH1D *>(inHisto.Get("h_all_cut"));
	auto h_all_cut_ce = static_cast<TH1D *>(inHisto.Get("h_all_cut_ce"));

	auto h_geometric_maxElayer_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_maxElayer_cut"));
	auto h_geometric_maxBarLayer_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_maxBarLayer_cut"));
	auto h_geometric_BGOTrackContainment_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_BGOTrackContainment_cut"));
	auto h_geometric_BGO_fiducial_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_BGO_fiducial_cut"));
	auto h_geometric_all_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_all_cut"));
	auto h_geometric_all_cut_ce = static_cast<TH1D *>(inHisto.Get("h_geometric_all_cut_ce"));

	auto h_BGOfiducial_nBarLayer13_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_nBarLayer13_cut"));
	auto h_BGOfiducial_maxRms_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_maxRms_cut"));
	auto h_BGOfiducial_track_selection_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_track_selection_cut"));
	auto h_BGOfiducial_psd_stk_match_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_psd_stk_match_cut"));
	auto h_BGOfiducial_psd_charge_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_psd_charge_cut"));
	auto h_BGOfiducial_stk_charge_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_stk_charge_cut"));
	auto h_BGOfiducial_xtrl_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_xtrl_cut"));
	auto h_BGOfiducial_all_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_all_cut"));
	auto h_BGOfiducial_all_cut_ce = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_all_cut_ce"));

	auto h_geo_BGOrec_slopeX = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOrec_slopeX"));
	auto h_geo_BGOrec_slopeY = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOrec_slopeY"));
	auto h_geo_BGOrec_interceptX = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOrec_interceptX"));
	auto h_geo_BGOrec_interceptY = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOrec_interceptY"));
	auto h_geo_BGOreco_topMap = static_cast<TH2D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOreco_topMap"));
	auto h_geo_BGOreco_bottomMap = static_cast<TH2D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOreco_bottomMap"));

	auto h_BGOrec_energy = static_cast<TH1D *>(inHisto.Get("BGO_Energy/h_BGOrec_energy"));
	auto h_layer_max_energy_ratio = static_cast<TH1D *>(inHisto.Get("BGO_Energy/h_layer_max_energy_ratio"));

	std::vector<TH1D *> h_layer_energy_ratio(DAMPE_bgo_nLayers);
	for (auto idx = 0; idx < DAMPE_bgo_nLayers; ++idx)
	{
		std::string histoName = "BGO_Energy/h_layer_energy_ratio_";
		std::ostringstream lNumber;
		lNumber << idx;
		histoName += lNumber.str();
		h_layer_energy_ratio[idx] = static_cast<TH1D *>(inHisto.Get(histoName.c_str()));
	}

	std::vector<TH2D*> sumRms_cosine (logEBins.size() - 1);
	for (auto it = sumRms_cosine.begin(); it != sumRms_cosine.end(); ++it)
	{
		std::string histo_name = "BGO_Energy/sumRms_cosine_" + to_string(std::distance(sumRms_cosine.begin(), it));
		(*it) = static_cast<TH2D *>(inHisto.Get(histo_name.c_str()));
	}

	auto sumRms_cosine_20_100 = static_cast<TH2D*>(inHisto.Get("BGO_Energy/sumRms_cosine_20_100"));
	auto sumRms_cosine_100_250 = static_cast<TH2D*>(inHisto.Get("BGO_Energy/sumRms_cosine_100_250"));
	auto sumRms_cosine_250_500 = static_cast<TH2D*>(inHisto.Get("BGO_Energy/sumRms_cosine_250_500"));
	auto sumRms_cosine_500_1000 = static_cast<TH2D*>(inHisto.Get("BGO_Energy/sumRms_cosine_500_1000"));
	auto sumRms_cosine_1000_3000 = static_cast<TH2D*>(inHisto.Get("BGO_Energy/sumRms_cosine_1000_3000"));
	auto sumRms_cosine_3000_10000 = static_cast<TH2D*>(inHisto.Get("BGO_Energy/sumRms_cosine_3000_10000"));

	auto h_xtrl_energy_int = static_cast<TH1D *>(inHisto.Get("xtrl/h_xtrl_energy_int"));
	auto h_xtrl = static_cast<TH2D *>(inHisto.Get("xtrl/h_xtrl"));
	auto e_discrimination_last = static_cast<TH2D *>(inHisto.Get("xtrl/e_discrimination_last"));
	auto e_discrimination_last_20_100 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_last_20_100"));
	auto e_discrimination_last_100_250 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_last_100_250"));
	auto e_discrimination_last_250_500 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_last_250_500"));
	auto e_discrimination_last_500_1000 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_last_500_1000"));
	auto e_discrimination_last_1000_3000 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_last_1000_3000"));
	auto e_discrimination_last_3000_10000 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_last_3000_10000"));
	auto e_discrimination = static_cast<TH2D *>(inHisto.Get("xtrl/e_discrimination"));
	auto e_discrimination_20_100 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_20_100"));
	auto e_discrimination_100_250 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_100_250"));
	auto e_discrimination_250_500 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_250_500"));
	auto e_discrimination_500_1000 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_500_1000"));
	auto e_discrimination_1000_3000 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_1000_3000"));
	auto e_discrimination_3000_10000 = static_cast<TH2D*>(inHisto.Get("xtrl/e_discrimination_3000_10000"));

	std::vector<TH1D*> bin_xtrl(logEBins.size() - 1);
	for (auto it = bin_xtrl.begin(); it != bin_xtrl.end(); ++it)
	{
		std::string histo_name = "xtrl/h_xtrl_bin_" + to_string(std::distance(bin_xtrl.begin(), it));
		(*it) = static_cast<TH1D *>(inHisto.Get(histo_name.c_str()));
	}

	auto h_psd_chargeX = static_cast<TH1D *>(inHisto.Get("PSDcharge/h_psd_chargeX"));
	auto h_psd_chargeY = static_cast<TH1D *>(inHisto.Get("PSDcharge/h_psd_chargeY"));
	auto h_psd_charge = static_cast<TH1D *>(inHisto.Get("PSDcharge/h_psd_charge"));
	auto h_psd_charge2D = static_cast<TH1D *>(inHisto.Get("PSDcharge/h_psd_charge2D"));

	auto h_psd_selected_chargeX = static_cast<TH1D *>(inHisto.Get("PSDcharge/h_psd_selected_chargeX"));
	auto h_psd_selected_chargeY = static_cast<TH1D *>(inHisto.Get("PSDcharge/h_psd_selected_chargeY"));
	auto h_psd_selected_charge = static_cast<TH1D *>(inHisto.Get("PSDcharge/h_psd_selected_charge"));
	auto h_psd_selected_charge2D = static_cast<TH1D *>(inHisto.Get("PSDcharge/h_psd_selected_charge2D"));

	auto h_stk_chargeX = static_cast<TH1D *>(inHisto.Get("STKcharge/h_stk_chargeX"));
	auto h_stk_chargeY = static_cast<TH1D *>(inHisto.Get("STKcharge/h_stk_chargeY"));
	auto h_stk_charge = static_cast<TH1D *>(inHisto.Get("STKcharge/h_stk_charge"));
	auto h_stk_charge2D = static_cast<TH1D *>(inHisto.Get("STKcharge/h_stk_charge2D"));

	auto h_stk_selected_chargeX = static_cast<TH1D *>(inHisto.Get("STKcharge/h_stk_selected_chargeX"));
	auto h_stk_selected_chargeY = static_cast<TH1D *>(inHisto.Get("STKcharge/h_stk_selected_chargeY"));
	auto h_stk_selected_charge = static_cast<TH1D *>(inHisto.Get("STKcharge/h_stk_selected_charge"));
	auto h_stk_selected_charge2D = static_cast<TH1D *>(inHisto.Get("STKcharge/h_stk_selected_charge2D"));

	auto h_background_under_xtrl_cut = static_cast<TH1D *>(inHisto.Get("mc_ancillary/h_background_under_xtrl_cut"));
	auto h_background_over_xtrl_cut = static_cast<TH1D *>(inHisto.Get("mc_ancillary/h_background_over_xtrl_cut"));

	h_trigger->SetDirectory(0);
	h_geometric_cut->SetDirectory(0);
	h_maxElayer_cut->SetDirectory(0);
	h_maxBarLayer_cut->SetDirectory(0);
	h_BGOTrackContainment_cut->SetDirectory(0);
	h_BGO_fiducial_cut->SetDirectory(0);
	h_all_cut->SetDirectory(0);
	h_all_cut_ce->SetDirectory(0);

	h_geometric_maxElayer_cut->SetDirectory(0);
	h_geometric_maxBarLayer_cut->SetDirectory(0);
	h_geometric_BGOTrackContainment_cut->SetDirectory(0);
	h_geometric_BGO_fiducial_cut->SetDirectory(0);
	h_geometric_all_cut->SetDirectory(0);
	h_geometric_all_cut_ce->SetDirectory(0);

	h_BGOfiducial_nBarLayer13_cut->SetDirectory(0);
	h_BGOfiducial_maxRms_cut->SetDirectory(0);
	h_BGOfiducial_track_selection_cut->SetDirectory(0);
	h_BGOfiducial_psd_stk_match_cut->SetDirectory(0);
	h_BGOfiducial_psd_charge_cut->SetDirectory(0);
	h_BGOfiducial_stk_charge_cut->SetDirectory(0);
	h_BGOfiducial_xtrl_cut->SetDirectory(0);
	h_BGOfiducial_all_cut->SetDirectory(0);
	h_BGOfiducial_all_cut_ce->SetDirectory(0);
	
	h_geo_BGOrec_slopeX->SetDirectory(0);
	h_geo_BGOrec_slopeY->SetDirectory(0);
	h_geo_BGOrec_interceptX->SetDirectory(0);
	h_geo_BGOrec_interceptY->SetDirectory(0);
	h_geo_BGOreco_topMap->SetDirectory(0);
	h_geo_BGOreco_bottomMap->SetDirectory(0);

	h_BGOrec_energy->SetDirectory(0);
	h_layer_max_energy_ratio->SetDirectory(0);

	for (auto idx = 0; idx < DAMPE_bgo_nLayers; ++idx)
		h_layer_energy_ratio[idx]->SetDirectory(0);

	for (auto it = sumRms_cosine.begin(); it != sumRms_cosine.end(); ++it)
		(*it)->SetDirectory(0);
	
	sumRms_cosine_20_100->SetDirectory(0);
	sumRms_cosine_100_250->SetDirectory(0);
	sumRms_cosine_250_500->SetDirectory(0);
	sumRms_cosine_500_1000->SetDirectory(0);
	sumRms_cosine_1000_3000->SetDirectory(0);
	sumRms_cosine_3000_10000->SetDirectory(0);

	h_xtrl_energy_int->SetDirectory(0);
	h_xtrl->SetDirectory(0);
	for (auto it = bin_xtrl.begin(); it != bin_xtrl.end(); ++it)
		(*it)->SetDirectory(0);

	e_discrimination_last->SetDirectory(0);
	e_discrimination_last_20_100->SetDirectory(0);
	e_discrimination_last_100_250->SetDirectory(0);
	e_discrimination_last_250_500->SetDirectory(0);
	e_discrimination_last_500_1000->SetDirectory(0);
	e_discrimination_last_1000_3000->SetDirectory(0);
	e_discrimination_last_3000_10000->SetDirectory(0);
	e_discrimination->SetDirectory(0);
	e_discrimination_20_100->SetDirectory(0);
	e_discrimination_100_250->SetDirectory(0);
	e_discrimination_250_500->SetDirectory(0);
	e_discrimination_500_1000->SetDirectory(0);
	e_discrimination_1000_3000->SetDirectory(0);
	e_discrimination_3000_10000->SetDirectory(0);

	h_psd_chargeX->SetDirectory(0);
	h_psd_chargeY->SetDirectory(0);
	h_psd_charge->SetDirectory(0);
	h_psd_charge2D->SetDirectory(0);

	h_psd_selected_chargeX->SetDirectory(0);
	h_psd_selected_chargeY->SetDirectory(0);
	h_psd_selected_charge->SetDirectory(0);
	h_psd_selected_charge2D->SetDirectory(0);

	h_stk_chargeX->SetDirectory(0);
	h_stk_chargeY->SetDirectory(0);
	h_stk_charge->SetDirectory(0);
	h_stk_charge2D->SetDirectory(0);

	h_stk_selected_chargeX->SetDirectory(0);
	h_stk_selected_chargeY->SetDirectory(0);
	h_stk_selected_charge->SetDirectory(0);
	h_stk_selected_charge2D->SetDirectory(0);

	h_background_under_xtrl_cut->SetDirectory(0);
	h_background_over_xtrl_cut->SetDirectory(0);

	inHisto.Close();

	// Write output TFile
	TFile outFile(outputPath.c_str(), "RECREATE");
	if (!outFile.IsOpen())
	{
		std::cerr << "\n\nError writing output TFile: " << outputPath << std::endl;
		exit(123);
	}

	// Write histos to file
	// Cut histos
	h_trigger->Write();
	h_geometric_cut->Write();
	h_maxElayer_cut->Write();
	h_maxBarLayer_cut->Write();
	h_BGOTrackContainment_cut->Write();
	h_BGO_fiducial_cut->Write();
	h_all_cut->Write();
	h_all_cut_ce->Write();

	// Cuts && Geometric Cut
	h_geometric_maxElayer_cut->Write();
	h_geometric_maxBarLayer_cut->Write();
	h_geometric_BGOTrackContainment_cut->Write();
	h_geometric_BGO_fiducial_cut->Write();
	h_geometric_all_cut->Write();
	h_geometric_all_cut_ce->Write();
	
	// Cuts && BGO fiducial volume cut
	h_BGOfiducial_nBarLayer13_cut->Write();
	h_BGOfiducial_maxRms_cut->Write();
	h_BGOfiducial_track_selection_cut->Write();
	h_BGOfiducial_psd_stk_match_cut->Write();
	h_BGOfiducial_psd_charge_cut->Write();
	h_BGOfiducial_stk_charge_cut->Write();
	h_BGOfiducial_xtrl_cut->Write();
	h_BGOfiducial_all_cut->Write();
	h_BGOfiducial_all_cut_ce->Write();

	// Create output ratio dir in the output TFile
	auto ratioDir = outFile.mkdir("Efficiency");

	// Create trigger folder
	auto trigger_dir = ratioDir->mkdir("Trigger");
	trigger_dir->cd();

	// Define TEfficiency pointers
	std::shared_ptr<TEfficiency> tr_eff_gometric_cut;
	std::shared_ptr<TEfficiency> tr_eff_maxElayer_cut;
	std::shared_ptr<TEfficiency> tr_eff_maxBarLayer_cut;
	std::shared_ptr<TEfficiency> tr_eff_BGOTrackContainment_cut;
	std::shared_ptr<TEfficiency> tr_eff_BGO_fiducial_cut;
	std::shared_ptr<TEfficiency> tr_eff_all_cut;

	if (TEfficiency::CheckConsistency(*h_geometric_cut, *h_trigger))
		tr_eff_gometric_cut = std::make_shared<TEfficiency>(*h_geometric_cut, *h_trigger);

	if (TEfficiency::CheckConsistency(*h_maxElayer_cut, *h_trigger))
		tr_eff_maxElayer_cut = std::make_shared<TEfficiency>(*h_maxElayer_cut, *h_trigger);

	if (TEfficiency::CheckConsistency(*h_maxBarLayer_cut, *h_trigger))
		tr_eff_maxBarLayer_cut = std::make_shared<TEfficiency>(*h_maxBarLayer_cut, *h_trigger);

	if (TEfficiency::CheckConsistency(*h_BGOTrackContainment_cut, *h_trigger))
		tr_eff_BGOTrackContainment_cut = std::make_shared<TEfficiency>(*h_BGOTrackContainment_cut, *h_trigger);

	if (TEfficiency::CheckConsistency(*h_BGO_fiducial_cut, *h_trigger))
		tr_eff_BGO_fiducial_cut = std::make_shared<TEfficiency>(*h_BGO_fiducial_cut, *h_trigger);

	if (TEfficiency::CheckConsistency(*h_all_cut, *h_trigger))
		tr_eff_all_cut = std::make_shared<TEfficiency>(*h_all_cut, *h_trigger);

	// Set uniform statistic option
	tr_eff_gometric_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_maxElayer_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_maxBarLayer_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_BGOTrackContainment_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_BGO_fiducial_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);
	
	tr_eff_gometric_cut->SetName("tr_eff_gometric_cut");
	tr_eff_maxElayer_cut->SetName("tr_eff_maxElayer_cut");
	tr_eff_maxBarLayer_cut->SetName("tr_eff_maxBarLayer_cut");
	tr_eff_BGOTrackContainment_cut->SetName("tr_eff_BGOTrackContainment_cut");
	tr_eff_BGO_fiducial_cut->SetName("tr_eff_BGO_fiducial_cut");
	tr_eff_all_cut->SetName("tr_eff_all_cut");
	
	tr_eff_gometric_cut->SetTitle("Gometric cut efficiency");
	tr_eff_maxElayer_cut->SetTitle("maxElayer cut efficiency");
	tr_eff_maxBarLayer_cut->SetTitle("maxBarLayer cut efficiency");
	tr_eff_BGOTrackContainment_cut->SetTitle("BGOTrackContainment cut efficiency");
	tr_eff_all_cut->SetTitle("all cut efficiency");

	// Write histos to disk
	tr_eff_gometric_cut->Write();
	tr_eff_maxElayer_cut->Write();
	tr_eff_maxBarLayer_cut->Write();
	tr_eff_BGOTrackContainment_cut->Write();
	tr_eff_BGO_fiducial_cut->Write();
	tr_eff_all_cut->Write();

	// Create geometric folder
	auto geometric_dir = ratioDir->mkdir("Geometric");
	geometric_dir->cd();

	// Define TEfficiency pointers
	std::shared_ptr<TEfficiency> geo_eff_maxElayer_cut;
	std::shared_ptr<TEfficiency> geo_eff_maxBarLayer_cut;
	std::shared_ptr<TEfficiency> geo_eff_BGOTrackContainment_cut;
	std::shared_ptr<TEfficiency> geo_eff_BGO_fiducial;
	std::shared_ptr<TEfficiency> geo_eff_all_cut;

	if (TEfficiency::CheckConsistency(*h_geometric_maxElayer_cut, *h_geometric_cut))
		geo_eff_maxElayer_cut = std::make_shared<TEfficiency>(*h_geometric_maxElayer_cut, *h_geometric_cut);

	if (TEfficiency::CheckConsistency(*h_geometric_maxBarLayer_cut, *h_geometric_cut))
		geo_eff_maxBarLayer_cut = std::make_shared<TEfficiency>(*h_geometric_maxBarLayer_cut, *h_geometric_cut);

	if (TEfficiency::CheckConsistency(*h_geometric_BGOTrackContainment_cut, *h_geometric_cut))
		geo_eff_BGOTrackContainment_cut = std::make_shared<TEfficiency>(*h_geometric_BGOTrackContainment_cut, *h_geometric_cut);

	if (TEfficiency::CheckConsistency(*h_geometric_BGO_fiducial_cut, *h_geometric_cut))
		geo_eff_BGO_fiducial = std::make_shared<TEfficiency>(*h_geometric_BGO_fiducial_cut, *h_geometric_cut);

	if (TEfficiency::CheckConsistency(*h_geometric_all_cut, *h_geometric_cut))
		geo_eff_all_cut = std::make_shared<TEfficiency>(*h_geometric_all_cut, *h_geometric_cut);

	// Set uniform statistic option
	geo_eff_maxElayer_cut->SetStatisticOption(TEfficiency::kBUniform);
	geo_eff_maxBarLayer_cut->SetStatisticOption(TEfficiency::kBUniform);
	geo_eff_BGOTrackContainment_cut->SetStatisticOption(TEfficiency::kBUniform);
	geo_eff_BGO_fiducial->SetStatisticOption(TEfficiency::kBUniform);
	geo_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);

	geo_eff_maxElayer_cut->SetName("geo_eff_maxElayer_cut");
	geo_eff_maxBarLayer_cut->SetName("geo_eff_maxBarLayer_cut");
	geo_eff_BGOTrackContainment_cut->SetName("geo_eff_BGOTrackContainment_cut");
	geo_eff_BGO_fiducial->SetName("geo_eff_BGO_fiducial");
	geo_eff_all_cut->SetName("geo_eff_all_cut");

	geo_eff_maxElayer_cut->SetTitle("geometic maxElayer cut efficiency");
	geo_eff_maxBarLayer_cut->SetTitle("geometric maxBarLayer cut efficiency");
	geo_eff_BGOTrackContainment_cut->SetTitle("geometric BGOTrackContainment cut efficiency");
	geo_eff_BGO_fiducial->SetTitle("geometric BGO fiducial cut efficiency");
	geo_eff_all_cut->SetTitle("geometric all cut efficiency");

	//Write histos to disk
	geo_eff_maxElayer_cut->Write();
	geo_eff_maxBarLayer_cut->Write();
	geo_eff_BGOTrackContainment_cut->Write();
	geo_eff_BGO_fiducial->Write();
	geo_eff_all_cut->Write();

	
	// Create BGO_fiducial_volume folder
	auto BGOfiducial_dir = ratioDir->mkdir("BGO_fiducial_volume");
	BGOfiducial_dir->cd();

	// Define TEfficiency pointers
	std::shared_ptr<TEfficiency> BGOfiducial_eff_nBarLayer13_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_maxRms_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_track_selection_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_psd_stk_match_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_psd_charge_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_stk_charge_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_xtrl_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_all_cut;

	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_maxRms_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_rms_track_selection_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut;


	if (TEfficiency::CheckConsistency(*h_BGOfiducial_nBarLayer13_cut, *h_BGO_fiducial_cut))
		BGOfiducial_eff_nBarLayer13_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_nBarLayer13_cut, *h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(*h_BGOfiducial_maxRms_cut, *h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_maxRms_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_maxRms_cut, *h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(*h_BGOfiducial_track_selection_cut, *h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_rms_track_selection_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_track_selection_cut, *h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(*h_BGOfiducial_psd_stk_match_cut, *h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_psd_stk_match_cut, *h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(*h_BGOfiducial_psd_charge_cut, *h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_psd_charge_cut, *h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(*h_BGOfiducial_stk_charge_cut, *h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_stk_charge_cut, *h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(*h_BGOfiducial_xtrl_cut, *h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_xtrl_cut, *h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(*h_BGOfiducial_all_cut, *h_BGO_fiducial_cut))
		BGOfiducial_eff_all_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_all_cut, *h_BGO_fiducial_cut);

	if (active_cuts.maxRms && active_cuts.nBarLayer13)
		if (TEfficiency::CheckConsistency(*h_BGOfiducial_maxRms_cut, *h_BGOfiducial_nBarLayer13_cut))
			BGOfiducial_eff_maxRms_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_maxRms_cut, *h_BGOfiducial_nBarLayer13_cut);

	if (active_cuts.track_selection && active_cuts.maxRms)
		if (TEfficiency::CheckConsistency(*h_BGOfiducial_track_selection_cut, *h_BGOfiducial_maxRms_cut))
			BGOfiducial_eff_track_selection_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_track_selection_cut, *h_BGOfiducial_maxRms_cut);

	if (active_cuts.psd_stk_match && active_cuts.track_selection)
		if (TEfficiency::CheckConsistency(*h_BGOfiducial_psd_stk_match_cut, *h_BGOfiducial_track_selection_cut))
			BGOfiducial_eff_psd_stk_match_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_psd_stk_match_cut, *h_BGOfiducial_track_selection_cut);

	if (active_cuts.psd_charge && active_cuts.psd_stk_match)
		if (TEfficiency::CheckConsistency(*h_BGOfiducial_psd_charge_cut, *h_BGOfiducial_psd_stk_match_cut))
			BGOfiducial_eff_psd_charge_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_psd_charge_cut, *h_BGOfiducial_psd_stk_match_cut);

	if (active_cuts.stk_charge && active_cuts.psd_charge)
		if (TEfficiency::CheckConsistency(*h_BGOfiducial_stk_charge_cut, *h_BGOfiducial_psd_stk_match_cut))
			BGOfiducial_eff_stk_charge_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_stk_charge_cut, *h_BGOfiducial_psd_stk_match_cut);

	if (active_cuts.xtrl && active_cuts.stk_charge)
		if (TEfficiency::CheckConsistency(*h_BGOfiducial_xtrl_cut, *h_BGOfiducial_stk_charge_cut))
			BGOfiducial_eff_xtrl_cut = std::make_shared<TEfficiency>(*h_BGOfiducial_xtrl_cut, *h_BGOfiducial_stk_charge_cut);

	// Set uniform statistic option
	BGOfiducial_eff_l13_maxRms_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_l13_rms_track_selection_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut->SetStatisticOption(TEfficiency::kBUniform);

	BGOfiducial_eff_nBarLayer13_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_maxRms_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_track_selection_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_psd_stk_match_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_psd_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_stk_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_xtrl_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);

	BGOfiducial_eff_l13_maxRms_cut->SetName("BGOfiducial_eff_l13_maxRms_cut");
	BGOfiducial_eff_l13_rms_track_selection_cut->SetName("BGOfiducial_eff_l13_rms_track_selection_cut");
	BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut->SetName("BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut");
	BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut->SetName("BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut");
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut->SetName("BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut");
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut->SetName("BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut");

	BGOfiducial_eff_nBarLayer13_cut->SetName("BGOfiducial_eff_nBarLayer13_cut");
	BGOfiducial_eff_maxRms_cut->SetName("BGOfiducial_eff_maxRms_cut");
	BGOfiducial_eff_track_selection_cut->SetName("BGOfiducial_eff_track_selection_cut");
	BGOfiducial_eff_psd_stk_match_cut->SetName("BGOfiducial_eff_psd_stk_match_cut");
	BGOfiducial_eff_psd_charge_cut->SetName("BGOfiducial_eff_psd_charge_cut");
	BGOfiducial_eff_stk_charge_cut->SetName("BGOfiducial_eff_stk_charge_cut");
	BGOfiducial_eff_xtrl_cut->SetName("BGOfiducial_eff_xtrl_cut");
	BGOfiducial_eff_all_cut->SetName("BGOfiducial_eff_all_cut");

	BGOfiducial_eff_l13_maxRms_cut->SetTitle("BGOfiducial l13 + maxRms cut efficiency");
	BGOfiducial_eff_l13_rms_track_selection_cut->SetTitle("BGOfiducial l13 + rms + track_selection cut efficiency");
	BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut->SetTitle("BGOfiducial l13 + rms + ts + PSD-STK match cut efficiency");
	BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut->SetTitle("BGOfiducial l13 + rms + ts + psdstk + PSD charge cut efficiency");
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut->SetTitle("BGOfiducial l13 + rms + ts + psdstk + pc + STK charge cut efficiency");
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut->SetTitle("BGOfiducial l13 + rms + ts + psdstk + pc + sc + xtrl cut efficiency");

	BGOfiducial_eff_nBarLayer13_cut->SetTitle("BGOfiducial nBarLayer13 cut efficiency");
	BGOfiducial_eff_maxRms_cut->SetTitle("BGOfiducial maxRms cut efficiency");
	BGOfiducial_eff_track_selection_cut->SetTitle("BGOfiducial track selection cut efficiency");
	BGOfiducial_eff_psd_stk_match_cut->SetTitle("BGOfiducial PSD-STK match cut efficiency");
	BGOfiducial_eff_psd_charge_cut->SetTitle("BGOfiducial PSD charge cut efficiency");
	BGOfiducial_eff_stk_charge_cut->SetTitle("BGOfiducial STK charge cut efficiency");
	BGOfiducial_eff_xtrl_cut->SetTitle("BGOfiducial xtrl cut efficiency");
	BGOfiducial_eff_all_cut->SetTitle("BGOfiducial all cut efficiency");

	// Write histos to disk
	BGOfiducial_eff_l13_maxRms_cut->Write();
	BGOfiducial_eff_l13_rms_track_selection_cut->Write();
	BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut->Write();
	BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut->Write();
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut->Write();
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut->Write();

	BGOfiducial_eff_nBarLayer13_cut->Write();
	BGOfiducial_eff_maxRms_cut->Write();
	BGOfiducial_eff_track_selection_cut->Write();
	BGOfiducial_eff_psd_stk_match_cut->Write();
	BGOfiducial_eff_psd_charge_cut->Write();
	BGOfiducial_eff_stk_charge_cut->Write();
	BGOfiducial_eff_xtrl_cut->Write();
	BGOfiducial_eff_all_cut->Write();
	
	auto geo_analysisDir = outFile.mkdir("Analysis_GeoCut");
	geo_analysisDir->cd();

	h_geo_BGOrec_slopeX->Write();
	h_geo_BGOrec_slopeY->Write();
	h_geo_BGOrec_interceptX->Write();
	h_geo_BGOrec_interceptY->Write();
	h_geo_BGOreco_topMap->Write();
	h_geo_BGOreco_bottomMap->Write();
	
	auto BGOdir = outFile.mkdir("BGO_Energy");
	BGOdir->cd();

	h_BGOrec_energy->Write();
	h_layer_max_energy_ratio->Write();

	for (auto it = sumRms_cosine.begin(); it != sumRms_cosine.end(); ++it)
		(*it)->Write();

	sumRms_cosine_20_100->Write();
	sumRms_cosine_100_250->Write();
	sumRms_cosine_250_500->Write();
	sumRms_cosine_500_1000->Write();
	sumRms_cosine_1000_3000->Write();
	sumRms_cosine_3000_10000->Write();

	for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
		h_layer_energy_ratio[lIdx]->Write();

	auto XTRLdir = outFile.mkdir("xtrl");
	XTRLdir->cd();

	h_xtrl_energy_int->Write();
	h_xtrl->Write();
	e_discrimination->Write();
	e_discrimination_20_100->Write();
	e_discrimination_100_250->Write();
	e_discrimination_250_500->Write();
	e_discrimination_500_1000->Write();
	e_discrimination_1000_3000->Write();
	e_discrimination_3000_10000->Write();
	e_discrimination_last->Write();
	e_discrimination_last_20_100->Write();
	e_discrimination_last_100_250->Write();
	e_discrimination_last_250_500->Write();
	e_discrimination_last_500_1000->Write();
	e_discrimination_last_1000_3000->Write();
	e_discrimination_last_3000_10000->Write();

	for (auto it = bin_xtrl.begin(); it != bin_xtrl.end(); ++it)
		(*it)->Write();

	auto psd_chargeDir = outFile.mkdir("PSDcharge");
	psd_chargeDir->cd();

	h_psd_chargeX->Write();
	h_psd_chargeY->Write();
	h_psd_charge->Write();
	h_psd_charge2D->Write();

	h_psd_selected_chargeX->Write();
	h_psd_selected_chargeY->Write();
	h_psd_selected_charge->Write();
	h_psd_selected_charge2D->Write();

	auto chargeDir = outFile.mkdir("STKcharge");
	chargeDir->cd();

	h_stk_chargeX->Write();
	h_stk_chargeY->Write();
	h_stk_charge->Write();
	h_stk_charge2D->Write();

	h_stk_selected_chargeX->Write();
	h_stk_selected_chargeY->Write();
	h_stk_selected_charge->Write();
	h_stk_selected_charge2D->Write();

	auto ancillaryDir = outFile.mkdir("proton_background");
	ancillaryDir->cd();

	h_background_under_xtrl_cut->Write();
	h_background_over_xtrl_cut->Write();

	// Create proton background ratio
	auto proton_background_ratio = static_cast<TH1D *>(h_background_under_xtrl_cut->Clone("proton_background_ratio"));
	proton_background_ratio->SetTitle("Proton background ratio");
	proton_background_ratio->Divide(h_background_over_xtrl_cut);

	proton_background_ratio->Write();

	// Close output TFile
	outFile.Close();
}