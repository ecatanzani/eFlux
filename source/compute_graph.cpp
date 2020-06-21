#include "binning.h"
#include "acceptance.h"
#include "wtsydp.h"
#include "read_sets_config_file.h"

#include "TH2D.h"
#include "TGraphErrors.h"
#include "TEfficiency.h"

void generateFinalGraph(
    const bool verbose,
    const bool pedantic,
    const std::string outputPath,
    const std::string complete_histo_path,
    const std::string wd)
{
    // Read energy log-binning
    auto logEBins = readLogBinning(wd);
    if (pedantic)
    {
        std::cout << "\nEnergy log binning..." << std::scientific;
        for (auto it = logEBins.begin(); it != logEBins.end(); ++it)
            std::cout << "\n"
                      << *it;
        std::cout << std::defaultfloat;
    }

    // Create and load acceptance events cuts from config file
    cuts_conf acceptance_cuts;
    data_active_cuts active_cuts;
    load_acceptance_struct(acceptance_cuts, active_cuts, wd);

    // Read dataSets connfig file
    data_set_conf input_sets;
    load_input_dsets_config(input_sets, wd);

    // Open input complete histos
    TFile inHisto(complete_histo_path.c_str(), "READ");
    if (!inHisto.IsOpen())
    {
        std::cerr << "Error reading input complete histo: " << complete_histo_path;
        exit(123);
    }
    
    auto h_geo_factor = static_cast<TH1D *>(inHisto.Get("h_geo_factor"));
    auto h_incoming = static_cast<TH1D *>(inHisto.Get("h_incoming"));
    auto h_trigger = static_cast<TH1D *>(inHisto.Get("h_trigger"));
    auto h_gometric_cut = static_cast<TH1D *>(inHisto.Get("h_gometric_cut"));
    auto h_maxElayer_cut = static_cast<TH1D *>(inHisto.Get("h_maxElayer_cut"));
    auto h_maxBarLayer_cut = static_cast<TH1D *>(inHisto.Get("h_maxBarLayer_cut"));
    auto h_BGOTrackContainment_cut = static_cast<TH1D *>(inHisto.Get("h_BGOTrackContainment_cut"));
    auto h_BGO_fiducial_cut = static_cast<TH1D *>(inHisto.Get("h_BGO_fiducial_cut"));
    auto h_nBarLayer13_cut = static_cast<TH1D *>(inHisto.Get("h_nBarLayer13_cut"));
    auto h_maxRms_cut = static_cast<TH1D *>(inHisto.Get("h_maxRms_cut"));
    auto h_track_selection_cut = static_cast<TH1D *>(inHisto.Get("h_track_selection_cut"));
    auto h_xtrl_cut = static_cast<TH1D *>(inHisto.Get("h_xtrl_cut"));
    auto h_psd_charge_cut = static_cast<TH1D *>(inHisto.Get("h_psd_charge_cut"));
    auto h_all_cut = static_cast<TH1D *>(inHisto.Get("h_all_cut"));

    auto h_geometric_maxElayer_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_maxElayer_cut"));
    auto h_geometric_maxBarLayer_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_maxBarLayer_cut"));
    auto h_geometric_BGOTrackContainment_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_BGOTrackContainment_cut"));
    auto h_geometric_BGO_fiducial_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_BGO_fiducial_cut"));
    auto h_geometric_nBarLayer13_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_nBarLayer13_cut"));
    auto h_geometric_maxRms_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_maxRms_cut"));
    auto h_geometric_track_selection_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_track_selection_cut"));
    auto h_geometric_xtrl_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_xtrl_cut"));
    auto h_geometric_psd_charge_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_psd_charge_cut"));
    auto h_geometric_all_cut = static_cast<TH1D *>(inHisto.Get("h_geometric_all_cut"));

    auto h_BGOfiducial_nBarLayer13_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_nBarLayer13_cut"));
    auto h_BGOfiducial_maxRms_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_maxRms_cut"));
    auto h_BGOfiducial_track_selection_cut = static_cast<TH1D*>(inHisto.Get("h_BGOfiducial_track_selection_cut"));
    auto h_BGOfiducial_xtrl_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_xtrl_cut"));
    auto h_BGOfiducial_psd_charge_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_psd_charge_cut"));
    auto h_BGOfiducial_all_cut = static_cast<TH1D *>(inHisto.Get("h_BGOfiducial_all_cut"));

    auto h_preGeo_BGOrec_topX_vs_realX = static_cast<TH1D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_BGOrec_topX_vs_realX"));
    auto h_preGeo_BGOrec_topY_vs_realY = static_cast<TH1D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_BGOrec_topY_vs_realY"));
    auto h_preGeo_real_slopeX = static_cast<TH1D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_real_slopeX"));
    auto h_preGeo_real_slopeY = static_cast<TH1D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_real_slopeY"));
    auto h_preGeo_BGOrec_slopeX = static_cast<TH1D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_BGOrec_slopeX"));
    auto h_preGeo_BGOrec_slopeY = static_cast<TH1D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_BGOrec_slopeY"));
    auto h_preGeo_real_interceptX = static_cast<TH1D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_real_interceptX"));
    auto h_preGeo_real_interceptY = static_cast<TH1D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_real_interceptY"));
    auto h_preGeo_BGOrec_interceptX = static_cast<TH1D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_BGOrec_interceptX"));
    auto h_preGeo_BGOrec_interceptY = static_cast<TH1D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_BGOrec_interceptY"));
    auto h_preGeo_real_topMap = static_cast<TH2D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_real_topMap"));
    auto h_preGeo_BGOreco_topMap = static_cast<TH2D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_BGOreco_topMap"));
    auto h_preGeo_real_bottomMap = static_cast<TH2D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_real_bottomMap"));
    auto h_preGeo_BGOreco_bottomMap = static_cast<TH2D *>(inHisto.Get("Analysis_preGeoCut/h_preGeo_BGOreco_bottomMap"));
    auto h_noBGOenergy_real_topMap = static_cast<TH2D *>(inHisto.Get("Analysis_preGeoCut/h_noBGOenergy_real_topMap"));

    auto h_geo_BGOrec_topX_vs_realX = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOrec_topX_vs_realX"));
    auto h_geo_BGOrec_topY_vs_realY = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOrec_topY_vs_realY"));
    auto h_geo_real_slopeX = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_real_slopeX"));
    auto h_geo_real_slopeY = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_real_slopeY"));
    auto h_geo_BGOrec_slopeX = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOrec_slopeX"));
    auto h_geo_BGOrec_slopeY = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOrec_slopeY"));
    auto h_geo_real_interceptX = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_real_interceptX"));
    auto h_geo_real_interceptY = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_real_interceptY"));
    auto h_geo_BGOrec_interceptX = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOrec_interceptX"));
    auto h_geo_BGOrec_interceptY = static_cast<TH1D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOrec_interceptY"));
    auto h_geo_real_topMap = static_cast<TH2D *>(inHisto.Get("Analysis_GeoCut/h_geo_real_topMap"));
    auto h_geo_BGOreco_topMap = static_cast<TH2D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOreco_topMap"));
    auto h_geo_real_bottomMap = static_cast<TH2D *>(inHisto.Get("Analysis_GeoCut/h_geo_real_bottomMap"));
    auto h_geo_BGOreco_bottomMap = static_cast<TH2D *>(inHisto.Get("Analysis_GeoCut/h_geo_BGOreco_bottomMap"));

    auto h_BGOrec_E = static_cast<TH1D *>(inHisto.Get("BGO_Energy/h_BGOrec_E"));
    auto h_BGOrec_E_corr = static_cast<TH1D *>(inHisto.Get("BGO_Energy/h_BGOrec_E_corr"));
    auto h_simu_energy = static_cast<TH1D *>(inHisto.Get("BGO_Energy/h_simu_energy"));
    auto h_energy_diff = static_cast<TH1D *>(inHisto.Get("BGO_Energy/h_energy_diff"));
    auto h_triggered_BGOrec_E = static_cast<TH1D *>(inHisto.Get("BGO_Energy/h_triggered_BGOrec_E"));
    auto h_triggered_BGOrec_E_corr = static_cast<TH1D *>(inHisto.Get("BGO_Energy/h_triggered_BGOrec_E_corr"));
    auto h_triggered_simu_energy = static_cast<TH1D *>(inHisto.Get("BGO_Energy/h_triggered_simu_energy"));
    auto h_triggered_energy_diff = static_cast<TH1D *>(inHisto.Get("BGO_Energy/h_triggered_energy_diff"));
    auto h_layer_max_energy_ratio = static_cast<TH1D *>(inHisto.Get("BGO_Energy/h_layer_max_energy_ratio"));
        
    std::vector<TH1D*> h_layer_energy_ratio(DAMPE_bgo_nLayers);
    for(auto idx=0; idx<DAMPE_bgo_nLayers; ++idx)
    {
        std::string histoName = "BGO_Energy/h_layer_energy_ratio_";
        std::ostringstream lNumber;
        lNumber << idx;
        histoName += lNumber.str();
        h_layer_energy_ratio[idx] = static_cast<TH1D *>(inHisto.Get(histoName.c_str()));
    }       

    auto h_xtrl_energy_int = static_cast<TH1D *>(inHisto.Get("xtrl/h_xtrl_energy_int"));
    auto h_xtrl = static_cast<TH2D *>(inHisto.Get("xtrl/h_xtrl"));

    h_geo_factor->SetDirectory(0);
    h_incoming->SetDirectory(0);
    h_trigger->SetDirectory(0);
    h_gometric_cut->SetDirectory(0);
    h_maxElayer_cut->SetDirectory(0);
    h_maxBarLayer_cut->SetDirectory(0);
    h_BGOTrackContainment_cut->SetDirectory(0);
    h_BGO_fiducial_cut->SetDirectory(0);
    h_nBarLayer13_cut->SetDirectory(0);
    h_maxRms_cut->SetDirectory(0);
    h_track_selection_cut->SetDirectory(0);
    h_xtrl_cut->SetDirectory(0);
    h_psd_charge_cut->SetDirectory(0);
    h_all_cut->SetDirectory(0);
    
    h_geometric_maxElayer_cut->SetDirectory(0);
    h_geometric_maxBarLayer_cut->SetDirectory(0);
    h_geometric_BGOTrackContainment_cut->SetDirectory(0);
    h_geometric_BGO_fiducial_cut->SetDirectory(0);
    h_geometric_nBarLayer13_cut->SetDirectory(0);
    h_geometric_maxRms_cut->SetDirectory(0);
    h_geometric_track_selection_cut->SetDirectory(0);
    h_geometric_xtrl_cut->SetDirectory(0);
    h_geometric_psd_charge_cut->SetDirectory(0);
    h_geometric_all_cut->SetDirectory(0);
    
    h_BGOfiducial_nBarLayer13_cut->SetDirectory(0);
    h_BGOfiducial_maxRms_cut->SetDirectory(0);
    h_BGOfiducial_track_selection_cut->SetDirectory(0);
    h_BGOfiducial_xtrl_cut->SetDirectory(0);
    h_BGOfiducial_psd_charge_cut->SetDirectory(0);
    h_BGOfiducial_all_cut->SetDirectory(0);
    
    h_preGeo_BGOrec_topX_vs_realX->SetDirectory(0);
    h_preGeo_BGOrec_topY_vs_realY->SetDirectory(0);
    h_preGeo_real_slopeX->SetDirectory(0);
    h_preGeo_real_slopeY->SetDirectory(0);
    h_preGeo_BGOrec_slopeX->SetDirectory(0);
    h_preGeo_BGOrec_slopeY->SetDirectory(0);
    h_preGeo_real_interceptX->SetDirectory(0);
    h_preGeo_real_interceptY->SetDirectory(0);
    h_preGeo_BGOrec_interceptX->SetDirectory(0);
    h_preGeo_BGOrec_interceptY->SetDirectory(0);
    h_preGeo_real_topMap->SetDirectory(0);
    h_preGeo_BGOreco_topMap->SetDirectory(0);
    h_preGeo_real_bottomMap->SetDirectory(0);
    h_preGeo_BGOreco_bottomMap->SetDirectory(0);
    h_noBGOenergy_real_topMap->SetDirectory(0);
    
    h_geo_BGOrec_topX_vs_realX->SetDirectory(0);
    h_geo_BGOrec_topY_vs_realY->SetDirectory(0);
    h_geo_real_slopeX->SetDirectory(0);
    h_geo_real_slopeY->SetDirectory(0);
    h_geo_BGOrec_slopeX->SetDirectory(0);
    h_geo_BGOrec_slopeY->SetDirectory(0);
    h_geo_real_interceptX->SetDirectory(0);
    h_geo_real_interceptY->SetDirectory(0);
    h_geo_BGOrec_interceptX->SetDirectory(0);
    h_geo_BGOrec_interceptY->SetDirectory(0);
    h_geo_real_topMap->SetDirectory(0);
    h_geo_BGOreco_topMap->SetDirectory(0);
    h_geo_real_bottomMap->SetDirectory(0);
    h_geo_BGOreco_bottomMap->SetDirectory(0);
    
    h_BGOrec_E->SetDirectory(0);
    h_BGOrec_E_corr->SetDirectory(0);
    h_simu_energy->SetDirectory(0);
    h_energy_diff->SetDirectory(0);
    h_triggered_BGOrec_E->SetDirectory(0);
    h_triggered_BGOrec_E_corr->SetDirectory(0);
    h_triggered_simu_energy->SetDirectory(0);
    h_triggered_energy_diff->SetDirectory(0);
    h_layer_max_energy_ratio->SetDirectory(0);

    for(auto idx=0; idx<DAMPE_bgo_nLayers; ++idx)
        h_layer_energy_ratio[idx]->SetDirectory(0);
        
    h_xtrl_energy_int->SetDirectory(0);
    h_xtrl->SetDirectory(0);

    inHisto.Close();

    double genSurface = 4 * TMath::Pi() * pow(acceptance_cuts.vertex_radius, 2) / 2;
    double scaleFactor = TMath::Pi() * genSurface;

    // Building acceptance histos
    auto h_acceptance_geometric_factor = static_cast<TH1D *>(h_geo_factor->Clone("h_acceptance_geometric_factor"));
    auto h_acceptance_gometric_cut = static_cast<TH1D *>(h_gometric_cut->Clone("h_acceptance_gometric_cut"));
    auto h_acceptance_maxElayer_cut = static_cast<TH1D *>(h_maxElayer_cut->Clone("h_acceptance_maxElayer_cut"));
    auto h_acceptance_maxBarLayer_cut = static_cast<TH1D *>(h_maxBarLayer_cut->Clone("h_acceptance_maxBarLayer_cut"));
    auto h_acceptance_BGOTrackContainment_cut = static_cast<TH1D *>(h_BGOTrackContainment_cut->Clone("h_acceptance_BGOTrackContainment_cut"));
    auto h_acceptance_BGO_fiducial_cut = static_cast<TH1D *>(h_BGO_fiducial_cut->Clone("h_acceptance_BGO_fiducial_cut"));
    auto h_acceptance_BGO_fiducial_nBarLayer13_cut = static_cast<TH1D *>(h_BGOfiducial_nBarLayer13_cut->Clone("h_acceptance_BGO_fiducial_nBarLayer13_cut"));
    auto h_acceptance_BGO_fiducial_maxRms_cut = static_cast<TH1D *>(h_BGOfiducial_maxRms_cut->Clone("h_acceptance_BGO_fiducial_maxRms_cut"));
    auto h_acceptance_BGO_fiducial_track_selection_cut = static_cast<TH1D *>(h_BGOfiducial_track_selection_cut->Clone("h_acceptance_BGO_fiducial_track_selection_cut"));
    auto h_acceptance_BGO_fiducial_xtrl_cut = static_cast<TH1D *>(h_BGOfiducial_xtrl_cut->Clone("h_acceptance_BGO_fiducial_xtrl_cut"));
    auto h_acceptance_BGO_fiducial_psd_charge_cut = static_cast<TH1D *>(h_BGOfiducial_psd_charge_cut->Clone("h_acceptance_BGO_fiducial_psd_charge_cut"));
    auto h_acceptance_nBarLayer13_cut = static_cast<TH1D *>(h_nBarLayer13_cut->Clone("h_acceptance_nBarLayer13_cut"));
    auto h_acceptance_maxRms_cut = static_cast<TH1D *>(h_maxRms_cut->Clone("h_acceptance_maxRms_cut"));
    auto h_acceptance_track_selection_cut = static_cast<TH1D *>(h_track_selection_cut->Clone("h_acceptance_track_selection_cut"));
    auto h_acceptance_xtrl_cut = static_cast<TH1D *>(h_xtrl_cut->Clone("h_acceptance_xtrl_cut"));
    auto h_acceptance_psd_charge_cut = static_cast<TH1D *>(h_psd_charge_cut->Clone("h_acceptance_psd_charge_cut"));
    auto h_acceptance_all_cut = static_cast<TH1D *>(h_all_cut->Clone("h_acceptance_all_cut"));

    h_acceptance_geometric_factor->Divide(h_incoming);
    h_acceptance_gometric_cut->Divide(h_incoming);
    h_acceptance_maxElayer_cut->Divide(h_incoming);
    h_acceptance_maxBarLayer_cut->Divide(h_incoming);
    h_acceptance_BGOTrackContainment_cut->Divide(h_incoming);
    h_acceptance_BGO_fiducial_cut->Divide(h_incoming);
    h_acceptance_BGO_fiducial_nBarLayer13_cut->Divide(h_incoming);
    h_acceptance_BGO_fiducial_maxRms_cut->Divide(h_incoming);
    h_acceptance_BGO_fiducial_track_selection_cut->Divide(h_incoming);
    h_acceptance_BGO_fiducial_xtrl_cut->Divide(h_incoming);
    h_acceptance_BGO_fiducial_psd_charge_cut->Divide(h_incoming);
    h_acceptance_nBarLayer13_cut->Divide(h_incoming);
    h_acceptance_maxRms_cut->Divide(h_incoming);
    h_acceptance_track_selection_cut->Divide(h_incoming);
    h_acceptance_xtrl_cut->Divide(h_incoming);
    h_acceptance_psd_charge_cut->Divide(h_incoming);
    h_acceptance_all_cut->Divide(h_incoming);

    h_acceptance_geometric_factor->Scale(scaleFactor);
    h_acceptance_gometric_cut->Scale(scaleFactor);
    h_acceptance_maxElayer_cut->Scale(scaleFactor);
    h_acceptance_maxBarLayer_cut->Scale(scaleFactor);
    h_acceptance_BGOTrackContainment_cut->Scale(scaleFactor);
    h_acceptance_BGO_fiducial_cut->Scale(scaleFactor);
    h_acceptance_BGO_fiducial_nBarLayer13_cut->Scale(scaleFactor);
    h_acceptance_BGO_fiducial_maxRms_cut->Scale(scaleFactor);
    h_acceptance_BGO_fiducial_track_selection_cut->Scale(scaleFactor);
    h_acceptance_BGO_fiducial_xtrl_cut->Scale(scaleFactor);
    h_acceptance_BGO_fiducial_psd_charge_cut->Scale(scaleFactor);
    h_acceptance_nBarLayer13_cut->Scale(scaleFactor);
    h_acceptance_maxRms_cut->Scale(scaleFactor);
    h_acceptance_track_selection_cut->Scale(scaleFactor);
    h_acceptance_xtrl_cut->Scale(scaleFactor);
    h_acceptance_psd_charge_cut->Scale(scaleFactor);
    h_acceptance_all_cut->Scale(scaleFactor);

    // Builing vectors
    std::vector<double> energyValues(h_incoming->GetXaxis()->GetNbins(), 0);

    std::vector<double> acceptanceValues_geometric_factor(energyValues.size(), 0);
    std::vector<double> acceptanceValues_gometric_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxElayer_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxBarLayer_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGOTrackContainment_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGO_fiducial_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGO_fiducial_nBarLayer13_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGO_fiducial_maxRms_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGO_fiducial_track_selection_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGO_fiducial_xtrl_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGO_fiducial_psd_charge_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_nBarLayer13_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxRms_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_track_selection_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_xtrl_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_psd_charge_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_all_cut(energyValues.size(), 0);

    //Building histo errors on energy and
    std::vector<double> acceptanceError_geometric_factor(h_incoming->GetXaxis()->GetNbins(), 0);
    std::vector<double> acceptanceError_gometric_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_maxElayer_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_maxBarLayer_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_BGOTrackContainment_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_BGO_fiducial_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_BGO_fiducial_nBarLayer13_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_BGO_fiducial_maxRms_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_BGO_fiducial_track_selection_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_BGO_fiducial_xtrl_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_BGO_fiducial_psd_charge_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_nBarLayer13_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_maxRms_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_track_selection_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_xtrl_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_psd_charge_cut(acceptanceError_geometric_factor.size(), 0);
    std::vector<double> acceptanceError_all_cut(acceptanceError_geometric_factor.size(), 0);

    std::vector<double> energyError(energyValues.size(), 0);

    for (auto it = logEBins.begin(); it != (logEBins.end() - 1); ++it)
    {
        auto index = std::distance(logEBins.begin(), it);
        energyValues[index] = wtsydp(*it, *(it + 1), getInputPowerLawIndex(*it, *(it + 1), input_sets));

        acceptanceValues_geometric_factor[index] = h_acceptance_geometric_factor->GetBinContent(index + 1);
        acceptanceValues_gometric_cut[index] = h_acceptance_gometric_cut->GetBinContent(index + 1);
        acceptanceValues_maxElayer_cut[index] = h_acceptance_maxElayer_cut->GetBinContent(index + 1);
        acceptanceValues_maxBarLayer_cut[index] = h_acceptance_maxBarLayer_cut->GetBinContent(index + 1);
        acceptanceValues_BGOTrackContainment_cut[index] = h_acceptance_BGOTrackContainment_cut->GetBinContent(index + 1);
        acceptanceValues_BGO_fiducial_cut[index] = h_acceptance_BGO_fiducial_cut->GetBinContent(index + 1);
        acceptanceValues_BGO_fiducial_nBarLayer13_cut[index] = h_acceptance_BGO_fiducial_nBarLayer13_cut->GetBinContent(index + 1);
        acceptanceValues_BGO_fiducial_maxRms_cut[index] = h_acceptance_BGO_fiducial_maxRms_cut->GetBinContent(index + 1);
        acceptanceValues_BGO_fiducial_track_selection_cut[index] = h_acceptance_BGO_fiducial_track_selection_cut->GetBinContent(index + 1);
        acceptanceValues_BGO_fiducial_xtrl_cut[index] = h_acceptance_BGO_fiducial_xtrl_cut->GetBinContent(index + 1);
        acceptanceValues_BGO_fiducial_psd_charge_cut[index] = h_acceptance_BGO_fiducial_psd_charge_cut->GetBinContent(index + 1);
        acceptanceValues_nBarLayer13_cut[index] = h_acceptance_nBarLayer13_cut->GetBinContent(index + 1);
        acceptanceValues_maxRms_cut[index] = h_acceptance_maxRms_cut->GetBinContent(index + 1);
        acceptanceValues_track_selection_cut[index] = h_acceptance_track_selection_cut->GetBinContent(index + 1);
        acceptanceValues_xtrl_cut[index] = h_acceptance_xtrl_cut->GetBinContent(index + 1);
        acceptanceValues_psd_charge_cut[index] = h_acceptance_psd_charge_cut->GetBinContent(index + 1);
        acceptanceValues_all_cut[index] = h_acceptance_all_cut->GetBinContent(index + 1);

        acceptanceError_geometric_factor[index] = h_acceptance_geometric_factor->GetBinError(index + 1);
        acceptanceError_gometric_cut[index] = h_acceptance_gometric_cut->GetBinError(index + 1);
        acceptanceError_maxElayer_cut[index] = h_acceptance_maxElayer_cut->GetBinError(index + 1);
        acceptanceError_maxBarLayer_cut[index] = h_acceptance_maxBarLayer_cut->GetBinError(index + 1);
        acceptanceError_BGOTrackContainment_cut[index] = h_acceptance_BGOTrackContainment_cut->GetBinError(index + 1);
        acceptanceError_BGO_fiducial_cut[index] = h_acceptance_BGO_fiducial_cut->GetBinError(index + 1);
        acceptanceError_BGO_fiducial_nBarLayer13_cut[index] = h_acceptance_BGO_fiducial_nBarLayer13_cut->GetBinError(index + 1);
        acceptanceError_BGO_fiducial_maxRms_cut[index] = h_acceptance_BGO_fiducial_maxRms_cut->GetBinError(index + 1);
        acceptanceError_BGO_fiducial_track_selection_cut[index] = h_acceptance_BGO_fiducial_track_selection_cut->GetBinError(index + 1);
        acceptanceError_BGO_fiducial_xtrl_cut[index] = h_acceptance_BGO_fiducial_xtrl_cut->GetBinError(index + 1);
        acceptanceError_BGO_fiducial_psd_charge_cut[index] = h_acceptance_BGO_fiducial_psd_charge_cut->GetBinError(index + 1);
        acceptanceError_nBarLayer13_cut[index] = h_acceptance_nBarLayer13_cut->GetBinError(index + 1);
        acceptanceError_maxRms_cut[index] = h_acceptance_maxRms_cut->GetBinError(index + 1);
        acceptanceError_track_selection_cut[index] = h_acceptance_track_selection_cut->GetBinError(index + 1);
        acceptanceError_xtrl_cut[index] = h_acceptance_xtrl_cut->GetBinError(index + 1);
        acceptanceError_psd_charge_cut[index] = h_acceptance_psd_charge_cut->GetBinError(index + 1);
        acceptanceError_all_cut[index] = h_acceptance_all_cut->GetBinError(index + 1);
    }

    // Building graphs
    TGraphErrors gr_acceptance_geometric_factor(energyValues.size(), &(energyValues[0]), &(acceptanceValues_geometric_factor[0]), &(energyError[0]), &(acceptanceError_geometric_factor[0]));
    TGraphErrors gr_acceptance_gometric_cut(energyValues.size(), &energyValues[0], &acceptanceValues_gometric_cut[0], &(energyError[0]), &(acceptanceError_gometric_cut[0]));
    TGraphErrors gr_acceptance_maxElayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxElayer_cut[0], &(energyError[0]), &(acceptanceError_maxElayer_cut[0]));
    TGraphErrors gr_acceptance_maxBarLayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxBarLayer_cut[0], &(energyError[0]), &(acceptanceError_maxBarLayer_cut[0]));
    TGraphErrors gr_acceptance_BGOTrackContainment_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGOTrackContainment_cut[0], &(energyError[0]), &(acceptanceError_BGOTrackContainment_cut[0]));
    TGraphErrors gr_acceptance_BGO_fiducial_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_cut[0]));
    TGraphErrors gr_acceptance_BGO_fiducial_nBarLayer13_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_nBarLayer13_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_nBarLayer13_cut[0]));
    TGraphErrors gr_acceptance_BGO_fiducial_maxRms_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_maxRms_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_maxRms_cut[0]));
    TGraphErrors gr_acceptance_BGO_fiducial_track_selection_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_track_selection_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_track_selection_cut[0]));
    TGraphErrors gr_acceptance_BGO_fiducial_xtrl_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_xtrl_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_xtrl_cut[0]));
    TGraphErrors gr_acceptance_BGO_fiducial_psd_charge_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_psd_charge_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_psd_charge_cut[0]));
    TGraphErrors gr_acceptance_nBarLayer13_cut(energyValues.size(), &energyValues[0], &acceptanceValues_nBarLayer13_cut[0], &(energyError[0]), &(acceptanceError_nBarLayer13_cut[0]));
    TGraphErrors gr_acceptance_maxRms_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxRms_cut[0], &(energyError[0]), &acceptanceError_maxRms_cut[0]);
    TGraphErrors gr_acceptance_track_selection_cut(energyValues.size(), &energyValues[0], &acceptanceValues_track_selection_cut[0], &(energyError[0]), &(acceptanceError_track_selection_cut[0]));
    TGraphErrors gr_acceptance_xtrl_cut(energyValues.size(), &energyValues[0], &acceptanceValues_xtrl_cut[0], &(energyError[0]),&(acceptanceError_xtrl_cut[0]));
    TGraphErrors gr_acceptance_psd_charge_cut(energyValues.size(), &energyValues[0], &acceptanceValues_psd_charge_cut[0], &(energyError[0]), &(acceptanceError_psd_charge_cut[0]));
    TGraphErrors gr_acceptance_all_cut(energyValues.size(), &energyValues[0], &acceptanceValues_all_cut[0], &(energyError[0]), &(acceptanceError_all_cut[0]));

    gr_acceptance_geometric_factor.SetName("gr_acceptance_geometric_factor");
    gr_acceptance_gometric_cut.SetName("gr_acceptance_gometric_cut");
    gr_acceptance_maxElayer_cut.SetName("gr_acceptance_maxElayer_cut");
    gr_acceptance_maxBarLayer_cut.SetName("gr_acceptance_maxBarLayer_cut");
    gr_acceptance_BGOTrackContainment_cut.SetName("gr_acceptance_BGOTrackContainment_cut");
    gr_acceptance_BGO_fiducial_cut.SetName("gr_acceptance_BGO_fiducial_cut");
    gr_acceptance_BGO_fiducial_nBarLayer13_cut.SetName("gr_acceptance_BGO_fiducial_nBarLayer13_cut");
    gr_acceptance_BGO_fiducial_maxRms_cut.SetName("gr_acceptance_BGO_fiducial_maxRms_cut");
    gr_acceptance_BGO_fiducial_track_selection_cut.SetName("gr_acceptance_BGO_fiducial_track_selection_cut");
    gr_acceptance_BGO_fiducial_xtrl_cut.SetName("gr_acceptance_BGO_fiducial_xtrl_cut");
    gr_acceptance_BGO_fiducial_psd_charge_cut.SetName("gr_acceptance_BGO_fiducial_psd_charge_cut");
    gr_acceptance_nBarLayer13_cut.SetName("gr_acceptance_nBarLayer13_cut");
    gr_acceptance_maxRms_cut.SetName("gr_acceptance_maxRms_cut");
    gr_acceptance_track_selection_cut.SetName("gr_acceptance_track_selection_cut");
    gr_acceptance_xtrl_cut.SetName("gr_acceptance_xtrl_cut");
    gr_acceptance_psd_charge_cut.SetName("gr_acceptance_psd_charge_cut");
    gr_acceptance_all_cut.SetName("gr_acceptance_all_cut");

    gr_acceptance_geometric_factor.SetTitle("Geometric Factor");
    gr_acceptance_gometric_cut.SetTitle("Acceptance - geometric cut");
    gr_acceptance_maxElayer_cut.SetTitle("Acceptance - maxElateral cut");
    gr_acceptance_maxBarLayer_cut.SetTitle("Acceptance - maxBarLayer cut");
    gr_acceptance_BGOTrackContainment_cut.SetTitle("Acceptance - BGOTrackContainment cut");
    gr_acceptance_BGO_fiducial_cut.SetTitle("Acceptance - BGO fiducial volume cut");
    gr_acceptance_BGO_fiducial_nBarLayer13_cut.SetTitle("Acceptance - BGO fiducial volume + nBarLayer13 cut");
    gr_acceptance_BGO_fiducial_maxRms_cut.SetTitle("Acceptance - BGO fiducial volume + maxRms cut");
    gr_acceptance_BGO_fiducial_track_selection_cut.SetTitle("Acceptance - BGO fiducial volume + track selection cut");
    gr_acceptance_BGO_fiducial_xtrl_cut.SetTitle("Acceptance - BGO fiducial volume + XTRL cut");
    gr_acceptance_BGO_fiducial_psd_charge_cut.SetTitle("Acceptance - BGO fiducial volume + PSD charge selection cut");
    gr_acceptance_nBarLayer13_cut.SetTitle("Acceptance - nBarLayer13 cut");
    gr_acceptance_maxRms_cut.SetTitle("Acceptance - maxRms cut");
    gr_acceptance_track_selection_cut.SetTitle("Acceptance - track selection cut");
    gr_acceptance_xtrl_cut.SetTitle("Acceptance - XTRL cut");
    gr_acceptance_psd_charge_cut.SetTitle("Acceptance - PSD charge selection cut");
    gr_acceptance_all_cut.SetTitle("Acceptance - all cut");
    
    // Write output TFile
    TFile outFile(outputPath.c_str(), "RECREATE");
    if (!outFile.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outputPath << std::endl;
        exit(123);
    }
    
    // Write histos to file
    // Acceptance - First-Cut histos
    h_geo_factor->Write();
    h_incoming->Write();
    h_trigger->Write();
    h_gometric_cut->Write();
    h_maxElayer_cut->Write();
    h_maxBarLayer_cut->Write();
    h_BGOTrackContainment_cut->Write();
    h_BGO_fiducial_cut->Write();
    h_nBarLayer13_cut->Write();
    h_maxRms_cut->Write();
    h_track_selection_cut->Write();
    h_xtrl_cut->Write();
    h_psd_charge_cut->Write();
    h_all_cut->Write();

    // Acceptance - Cuts && Geometric Cut
    h_geometric_maxElayer_cut->Write();
    h_geometric_maxBarLayer_cut->Write();
    h_geometric_BGOTrackContainment_cut->Write();
    h_geometric_BGO_fiducial_cut->Write();
    h_geometric_nBarLayer13_cut->Write();
    h_geometric_maxRms_cut->Write();
    h_geometric_track_selection_cut->Write();
    h_geometric_xtrl_cut->Write();
    h_geometric_psd_charge_cut->Write();
    h_geometric_all_cut->Write();
    // Acceptance - Cuts && BGO fiducial volume cut
    h_BGOfiducial_nBarLayer13_cut->Write();
    h_BGOfiducial_maxRms_cut->Write();
    h_BGOfiducial_track_selection_cut->Write();
    h_BGOfiducial_xtrl_cut->Write();
    h_BGOfiducial_psd_charge_cut->Write();
    h_BGOfiducial_all_cut->Write();

    // Create output acceptance_histo dir in the output TFile
    auto acceptanceHistoDir = outFile.mkdir("Acceptance_histos");
    acceptanceHistoDir->cd();

    h_acceptance_geometric_factor->Write();
    h_acceptance_gometric_cut->Write();
    h_acceptance_maxElayer_cut->Write();
    h_acceptance_maxBarLayer_cut->Write();
    h_acceptance_BGOTrackContainment_cut->Write();
    h_acceptance_BGO_fiducial_cut->Write();
    h_acceptance_BGO_fiducial_nBarLayer13_cut->Write();
    h_acceptance_BGO_fiducial_maxRms_cut->Write();
    h_acceptance_BGO_fiducial_track_selection_cut->Write();
    h_acceptance_BGO_fiducial_xtrl_cut->Write();
    h_acceptance_BGO_fiducial_psd_charge_cut->Write();
    h_acceptance_nBarLayer13_cut->Write();
    h_acceptance_maxRms_cut->Write();
    h_acceptance_track_selection_cut->Write();
    h_acceptance_xtrl_cut->Write();
    h_acceptance_psd_charge_cut->Write();
    h_acceptance_all_cut->Write();

    outFile.cd();

    // Create output acceptance dir in the output TFile
    auto acceptanceDir = outFile.mkdir("Acceptance");
    acceptanceDir->cd();

    // Write final TGraphs
    gr_acceptance_gometric_cut.Write();
    gr_acceptance_maxElayer_cut.Write();
    gr_acceptance_maxBarLayer_cut.Write();
    gr_acceptance_BGOTrackContainment_cut.Write();
    gr_acceptance_BGO_fiducial_cut.Write();
    gr_acceptance_BGO_fiducial_nBarLayer13_cut.Write();
    gr_acceptance_BGO_fiducial_maxRms_cut.Write();
    gr_acceptance_BGO_fiducial_track_selection_cut.Write();
    gr_acceptance_BGO_fiducial_xtrl_cut.Write();
    gr_acceptance_BGO_fiducial_psd_charge_cut.Write();
    gr_acceptance_nBarLayer13_cut.Write();
    gr_acceptance_maxRms_cut.Write();
    gr_acceptance_track_selection_cut.Write();
    gr_acceptance_xtrl_cut.Write();
    gr_acceptance_psd_charge_cut.Write();
    gr_acceptance_all_cut.Write();

    // Return to main TFile directory
    outFile.cd();
    
    auto geoFactor = outFile.mkdir("GeometricFactor");
    geoFactor->cd();
    
    gr_acceptance_geometric_factor.Write();

    outFile.cd();

    // Create output ratio dir in the output TFile
    auto ratioDir = outFile.mkdir("Efficiency");

    // Create trigger folder
    auto trigger_dir = ratioDir->mkdir("Trigger");
    trigger_dir->cd();

    // Define TEfficiency pointers
    TEfficiency *trigger_efficiency = nullptr;
    TEfficiency *tr_eff_gometric_cut = nullptr;
    TEfficiency *tr_eff_maxElayer_cut = nullptr;
    TEfficiency *tr_eff_maxBarLayer_cut = nullptr;
    TEfficiency *tr_eff_BGOTrackContainment_cut = nullptr;
    TEfficiency *tr_eff_BGO_fiducial_cut = nullptr;
    TEfficiency *tr_eff_nBarLayer13_cut = nullptr;
    TEfficiency *tr_eff_maxRms_cut = nullptr;
    TEfficiency *tr_eff_track_selection_cut = nullptr;
    TEfficiency *tr_eff_xtrl_cut = nullptr;
    TEfficiency *tr_eff_psd_charge_cut = nullptr;
    TEfficiency *tr_eff_all_cut = nullptr;

    if (TEfficiency::CheckConsistency(*h_gometric_cut, *h_gometric_cut))
        trigger_efficiency = new TEfficiency(*h_gometric_cut, *h_gometric_cut);

    if (TEfficiency::CheckConsistency(*h_gometric_cut, *h_trigger))
        tr_eff_gometric_cut = new TEfficiency(*h_gometric_cut, *h_trigger);

    if (TEfficiency::CheckConsistency(*h_maxElayer_cut, *h_trigger))
        tr_eff_maxElayer_cut = new TEfficiency(*h_maxElayer_cut, *h_trigger);

    if (TEfficiency::CheckConsistency(*h_maxBarLayer_cut, *h_trigger))
        tr_eff_maxBarLayer_cut = new TEfficiency(*h_maxBarLayer_cut, *h_trigger);

    if (TEfficiency::CheckConsistency(*h_BGOTrackContainment_cut, *h_trigger))
        tr_eff_BGOTrackContainment_cut = new TEfficiency(*h_BGOTrackContainment_cut, *h_trigger);

    if (TEfficiency::CheckConsistency(*h_BGO_fiducial_cut, *h_trigger))
        tr_eff_BGO_fiducial_cut = new TEfficiency(*h_BGO_fiducial_cut, *h_trigger);

    if (TEfficiency::CheckConsistency(*h_nBarLayer13_cut, *h_trigger))
        tr_eff_nBarLayer13_cut = new TEfficiency(*h_nBarLayer13_cut, *h_trigger);

    if (TEfficiency::CheckConsistency(*h_maxRms_cut, *h_trigger))
        tr_eff_maxRms_cut = new TEfficiency(*h_maxRms_cut, *h_trigger);

    if (TEfficiency::CheckConsistency(*h_track_selection_cut, *h_trigger))
        tr_eff_track_selection_cut = new TEfficiency(*h_track_selection_cut, *h_trigger);

    if (TEfficiency::CheckConsistency(*h_xtrl_cut, *h_trigger))
        tr_eff_xtrl_cut = new TEfficiency(*h_xtrl_cut, *h_trigger);

    if (TEfficiency::CheckConsistency(*h_psd_charge_cut, *h_trigger))
        tr_eff_psd_charge_cut = new TEfficiency(*h_psd_charge_cut, *h_trigger);

    if (TEfficiency::CheckConsistency(*h_all_cut, *h_trigger))
        tr_eff_all_cut = new TEfficiency(*h_all_cut, *h_trigger);

    // Set uniform statistic option
    trigger_efficiency->SetStatisticOption(TEfficiency::kBUniform);
    tr_eff_gometric_cut->SetStatisticOption(TEfficiency::kBUniform);
    tr_eff_maxElayer_cut->SetStatisticOption(TEfficiency::kBUniform);
    tr_eff_maxBarLayer_cut->SetStatisticOption(TEfficiency::kBUniform);
    tr_eff_BGOTrackContainment_cut->SetStatisticOption(TEfficiency::kBUniform);
    tr_eff_BGO_fiducial_cut->SetStatisticOption(TEfficiency::kBUniform);
    tr_eff_nBarLayer13_cut->SetStatisticOption(TEfficiency::kBUniform);
    tr_eff_maxRms_cut->SetStatisticOption(TEfficiency::kBUniform);
    tr_eff_track_selection_cut->SetStatisticOption(TEfficiency::kBUniform);
    tr_eff_xtrl_cut->SetStatisticOption(TEfficiency::kBUniform);
    tr_eff_psd_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
    tr_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);

    trigger_efficiency->SetName("trigger_efficiency");
    tr_eff_gometric_cut->SetName("tr_eff_gometric_cut");
    tr_eff_maxElayer_cut->SetName("tr_eff_maxElayer_cut");
    tr_eff_maxBarLayer_cut->SetName("tr_eff_maxBarLayer_cut");
    tr_eff_BGOTrackContainment_cut->SetName("tr_eff_BGOTrackContainment_cut");
    tr_eff_BGO_fiducial_cut->SetName("tr_eff_BGO_fiducial_cut");
    tr_eff_nBarLayer13_cut->SetName("tr_eff_nBarLayer13_cut");
    tr_eff_maxRms_cut->SetName("tr_eff_maxRms_cut");
    tr_eff_track_selection_cut->SetName("tr_eff_track_selection_cut");
    tr_eff_xtrl_cut->SetName("tr_eff_xtrl_cut");
    tr_eff_psd_charge_cut->SetName("tr_eff_psd_charge_cut");
    tr_eff_all_cut->SetName("tr_eff_all_cut");

    trigger_efficiency->SetTitle("Trigger efficiency");
    tr_eff_gometric_cut->SetTitle("Gometric cut efficiency");
    tr_eff_maxElayer_cut->SetTitle("maxElayer cut efficiency");
    tr_eff_maxBarLayer_cut->SetTitle("maxBarLayer cut efficiency");
    tr_eff_BGOTrackContainment_cut->SetTitle("BGOTrackContainment cut efficiency");
    tr_eff_BGO_fiducial_cut->SetTitle("BGO fiducial cut efficiency");
    tr_eff_nBarLayer13_cut->SetTitle("nBarLayer13 cut efficiency");
    tr_eff_maxRms_cut->SetTitle("maxRms cut efficiency");
    tr_eff_track_selection_cut->SetTitle("track selection cut efficiency");
    tr_eff_xtrl_cut->SetTitle("xtrl cut efficiency");
    tr_eff_psd_charge_cut->SetTitle("psd charge cut efficiency");
    tr_eff_all_cut->SetTitle("all cut efficiency");

    // Write histos to disk
    trigger_efficiency->Write();
    tr_eff_gometric_cut->Write();
    tr_eff_maxElayer_cut->Write();
    tr_eff_maxBarLayer_cut->Write();
    tr_eff_BGOTrackContainment_cut->Write();
    tr_eff_BGO_fiducial_cut->Write();
    tr_eff_nBarLayer13_cut->Write();
    tr_eff_maxRms_cut->Write();
    tr_eff_track_selection_cut->Write();
    tr_eff_xtrl_cut->Write();
    tr_eff_psd_charge_cut->Write();
    tr_eff_all_cut->Write();

    // Clean memory
    trigger_efficiency->Delete();
    tr_eff_gometric_cut->Delete();
    tr_eff_maxElayer_cut->Delete();
    tr_eff_maxBarLayer_cut->Delete();
    tr_eff_BGOTrackContainment_cut->Delete();
    tr_eff_BGO_fiducial_cut->Delete();
    tr_eff_nBarLayer13_cut->Delete();
    tr_eff_maxRms_cut->Delete();
    tr_eff_track_selection_cut->Delete();
    tr_eff_xtrl_cut->Delete();
    tr_eff_psd_charge_cut->Delete();
    tr_eff_all_cut->Delete();

    // Create geometric folder
    auto geometric_dir = ratioDir->mkdir("Geometric");
    geometric_dir->cd();

    // Define TEfficiency pointers
    TEfficiency *geo_eff_maxElayer_cut = nullptr;
    TEfficiency *geo_eff_maxBarLayer_cut = nullptr;
    TEfficiency *geo_eff_BGOTrackContainment_cut = nullptr;
    TEfficiency *geo_eff_BGO_fiducial = nullptr;
    TEfficiency *geo_eff_nBarLayer13_cut = nullptr;
    TEfficiency *geo_eff_maxRms_cut = nullptr;
    TEfficiency *geo_eff_track_selection_cut = nullptr;
    TEfficiency *geo_eff_xtrl_cut = nullptr;
    TEfficiency *geo_eff_psd_charge_cut = nullptr;
    TEfficiency *geo_eff_all_cut = nullptr;

    if (TEfficiency::CheckConsistency(*h_geometric_maxElayer_cut, *h_gometric_cut))
        geo_eff_maxElayer_cut = new TEfficiency(*h_geometric_maxElayer_cut, *h_gometric_cut);
    
    if (TEfficiency::CheckConsistency(*h_geometric_maxBarLayer_cut, *h_gometric_cut))
        geo_eff_maxBarLayer_cut = new TEfficiency(*h_geometric_maxBarLayer_cut, *h_gometric_cut);

    if (TEfficiency::CheckConsistency(*h_geometric_BGOTrackContainment_cut, *h_gometric_cut))
        geo_eff_BGOTrackContainment_cut = new TEfficiency(*h_geometric_BGOTrackContainment_cut, *h_gometric_cut);
    
    if (TEfficiency::CheckConsistency(*h_geometric_BGO_fiducial_cut, *h_gometric_cut))
        geo_eff_BGO_fiducial = new TEfficiency(*h_geometric_BGO_fiducial_cut, *h_gometric_cut);

    if (TEfficiency::CheckConsistency(*h_geometric_nBarLayer13_cut, *h_gometric_cut))
        geo_eff_nBarLayer13_cut = new TEfficiency(*h_geometric_nBarLayer13_cut, *h_gometric_cut);

    if (TEfficiency::CheckConsistency(*h_geometric_maxRms_cut, *h_gometric_cut))
        geo_eff_maxRms_cut = new TEfficiency(*h_geometric_maxRms_cut, *h_gometric_cut);
    
    if (TEfficiency::CheckConsistency(*h_geometric_track_selection_cut, *h_gometric_cut))
        geo_eff_track_selection_cut = new TEfficiency(*h_geometric_track_selection_cut, *h_gometric_cut);
    
    if (TEfficiency::CheckConsistency(*h_geometric_xtrl_cut, *h_gometric_cut))
        geo_eff_xtrl_cut = new TEfficiency(*h_geometric_xtrl_cut, *h_gometric_cut);

    if (TEfficiency::CheckConsistency(*h_geometric_psd_charge_cut, *h_gometric_cut))
        geo_eff_psd_charge_cut = new TEfficiency(*h_geometric_psd_charge_cut, *h_gometric_cut);
    
    if (TEfficiency::CheckConsistency(*h_geometric_all_cut, *h_gometric_cut))
        geo_eff_all_cut = new TEfficiency(*h_geometric_all_cut, *h_gometric_cut);
    
    // Set uniform statistic option
    geo_eff_maxElayer_cut->SetStatisticOption(TEfficiency::kBUniform);
    geo_eff_maxBarLayer_cut->SetStatisticOption(TEfficiency::kBUniform);
    geo_eff_BGOTrackContainment_cut->SetStatisticOption(TEfficiency::kBUniform);
    geo_eff_BGO_fiducial->SetStatisticOption(TEfficiency::kBUniform);
    geo_eff_nBarLayer13_cut->SetStatisticOption(TEfficiency::kBUniform);
    geo_eff_maxRms_cut->SetStatisticOption(TEfficiency::kBUniform);
    geo_eff_track_selection_cut->SetStatisticOption(TEfficiency::kBUniform);
    geo_eff_xtrl_cut->SetStatisticOption(TEfficiency::kBUniform);
    geo_eff_psd_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
    geo_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);

    //Write histos to disk
    geo_eff_maxElayer_cut->Write();
    geo_eff_maxBarLayer_cut->Write();
    geo_eff_BGOTrackContainment_cut->Write();
    geo_eff_BGO_fiducial->Write();
    geo_eff_nBarLayer13_cut->Write();
    geo_eff_maxRms_cut->Write();
    geo_eff_track_selection_cut->Write();
    geo_eff_xtrl_cut->Write();
    geo_eff_psd_charge_cut->Write();
    geo_eff_all_cut->Write();

    // Clean memory
    geo_eff_maxElayer_cut->Delete();
    geo_eff_maxBarLayer_cut->Delete();
    geo_eff_BGOTrackContainment_cut->Delete();
    geo_eff_BGO_fiducial->Delete();
    geo_eff_nBarLayer13_cut->Delete();
    geo_eff_maxRms_cut->Delete();
    geo_eff_track_selection_cut->Delete();
    geo_eff_xtrl_cut->Delete();
    geo_eff_psd_charge_cut->Delete();
    geo_eff_all_cut->Delete();

    auto BGOfiducial_dir = ratioDir->mkdir("BGO_fiducial_volume");
    BGOfiducial_dir->cd();

    // Define TEfficiency pointers
    TEfficiency *BGOfiducial_eff_nBarLayer13_cut = nullptr;
    TEfficiency *BGOfiducial_eff_maxRms_cut = nullptr;
    TEfficiency *BGOfiducial_eff_track_selection_cut = nullptr;
    TEfficiency *BGOfiducial_eff_xtrl_cut = nullptr;
    TEfficiency *BGOfiducial_eff_psd_charge_cut = nullptr;
    TEfficiency *BGOfiducial_eff_all_cut = nullptr;

    if (TEfficiency::CheckConsistency(*h_BGOfiducial_nBarLayer13_cut, *h_BGO_fiducial_cut))
        BGOfiducial_eff_nBarLayer13_cut = new TEfficiency(*h_BGOfiducial_nBarLayer13_cut, *h_BGO_fiducial_cut);

    if (TEfficiency::CheckConsistency(*h_BGOfiducial_maxRms_cut, *h_BGO_fiducial_cut))
        BGOfiducial_eff_maxRms_cut = new TEfficiency(*h_BGOfiducial_maxRms_cut, *h_BGO_fiducial_cut);

    if (TEfficiency::CheckConsistency(*h_BGOfiducial_track_selection_cut, *h_BGO_fiducial_cut))
        BGOfiducial_eff_track_selection_cut = new TEfficiency(*h_BGOfiducial_track_selection_cut, *h_BGO_fiducial_cut);

    if (TEfficiency::CheckConsistency(*h_BGOfiducial_xtrl_cut, *h_BGO_fiducial_cut))
        BGOfiducial_eff_xtrl_cut = new TEfficiency(*h_BGOfiducial_xtrl_cut, *h_BGO_fiducial_cut);

    if (TEfficiency::CheckConsistency(*h_BGOfiducial_psd_charge_cut, *h_BGO_fiducial_cut))
        BGOfiducial_eff_psd_charge_cut = new TEfficiency(*h_BGOfiducial_psd_charge_cut, *h_BGO_fiducial_cut);

    if (TEfficiency::CheckConsistency(*h_BGOfiducial_all_cut, *h_BGO_fiducial_cut))
        BGOfiducial_eff_all_cut = new TEfficiency(*h_BGOfiducial_all_cut, *h_BGO_fiducial_cut);
    
    // Set uniform statistic option
    BGOfiducial_eff_nBarLayer13_cut->SetStatisticOption(TEfficiency::kBUniform);
    BGOfiducial_eff_maxRms_cut->SetStatisticOption(TEfficiency::kBUniform);
    BGOfiducial_eff_track_selection_cut->SetStatisticOption(TEfficiency::kBUniform);
    BGOfiducial_eff_xtrl_cut->SetStatisticOption(TEfficiency::kBUniform);
    BGOfiducial_eff_psd_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
    BGOfiducial_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);

    // Write histos to disk
    BGOfiducial_eff_nBarLayer13_cut->Write();
    BGOfiducial_eff_maxRms_cut->Write();
    BGOfiducial_eff_track_selection_cut->Write();
    BGOfiducial_eff_xtrl_cut->Write();
    BGOfiducial_eff_psd_charge_cut->Write();
    BGOfiducial_eff_all_cut->Write();

    // Clean memory
    BGOfiducial_eff_nBarLayer13_cut->Delete();
    BGOfiducial_eff_maxRms_cut->Delete();
    BGOfiducial_eff_track_selection_cut->Delete();
    BGOfiducial_eff_xtrl_cut->Delete();
    BGOfiducial_eff_psd_charge_cut->Delete();
    BGOfiducial_eff_all_cut->Delete();
    
    // Create output analysis dir in the output TFile
    auto preGeo_analysisDir = outFile.mkdir("Analysis_preGeoCut");
    preGeo_analysisDir->cd();

    h_preGeo_BGOrec_topX_vs_realX->Write();
    h_preGeo_BGOrec_topY_vs_realY->Write();
    h_preGeo_real_slopeX->Write();
    h_preGeo_real_slopeY->Write();
    h_preGeo_BGOrec_slopeX->Write();
    h_preGeo_BGOrec_slopeY->Write();
    h_preGeo_real_interceptX->Write();
    h_preGeo_real_interceptY->Write();
    h_preGeo_BGOrec_interceptX->Write();
    h_preGeo_BGOrec_interceptY->Write();
    h_preGeo_real_topMap->Write();
    h_preGeo_BGOreco_topMap->Write();
    h_preGeo_real_bottomMap->Write();
    h_preGeo_BGOreco_bottomMap->Write();

    h_noBGOenergy_real_topMap->Write();

    outFile.cd();

    auto geo_analysisDir = outFile.mkdir("Analysis_GeoCut");
    geo_analysisDir->cd();

    h_geo_BGOrec_topX_vs_realX->Write();
    h_geo_BGOrec_topY_vs_realY->Write();
    h_geo_real_slopeX->Write();
    h_geo_real_slopeY->Write();
    h_geo_BGOrec_slopeX->Write();
    h_geo_BGOrec_slopeY->Write();
    h_geo_real_interceptX->Write();
    h_geo_real_interceptY->Write();
    h_geo_BGOrec_interceptX->Write();
    h_geo_BGOrec_interceptY->Write();
    h_geo_real_topMap->Write();
    h_geo_BGOreco_topMap->Write();
    h_geo_real_bottomMap->Write();
    h_geo_BGOreco_bottomMap->Write();

    outFile.cd();
    
    auto BGOdir = outFile.mkdir("BGO_Energy");
    BGOdir->cd();

    h_BGOrec_E->Write();
    h_BGOrec_E_corr->Write();
    h_simu_energy->Write();
    h_energy_diff->Write();
    h_layer_max_energy_ratio->Write();

    h_triggered_BGOrec_E->Write();
    h_triggered_BGOrec_E_corr->Write();
    h_triggered_simu_energy->Write();
    h_triggered_energy_diff->Write();
    
    for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
        h_layer_energy_ratio[lIdx]->Write();

    auto XTRLdir = outFile.mkdir("xtrl");
    XTRLdir->cd();

    h_xtrl_energy_int->Write();
    h_xtrl->Write();

    // Close output TFile
    outFile.Close();
    
}