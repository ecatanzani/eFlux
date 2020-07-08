#include "acceptance.h"
#include "acceptance_cuts.h"
#include "data_cuts.h"
#include "energy_match.h"
#include "aggregate_events.h"
#include "wtsydp.h"
#include "BGO_energy_cuts.h"
#include "DAMPE_geo_structure.h"
#include "DmpBgoContainer.h"
#include "read_sets_config_file.h"
#include "binning.h"
#include "charge.h"
#include "mc_ancillary.h"

#include "TGraphErrors.h"
#include "TEfficiency.h"

inline void init_BGO_histos(std::vector<TH1D> &h_layer_energy_ratio)
{
    h_layer_energy_ratio.resize(DAMPE_bgo_nLayers);

    for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
    {
        TString h_ratio_name = "h_layer_energy_ratio_";
        TString h_ratio_title = "Energy Ratio - BGO layer ";

        h_ratio_name += lIdx;
        h_ratio_title += lIdx;

        h_layer_energy_ratio[lIdx] = TH1D(h_ratio_name.Data(), h_ratio_title.Data(), 100, 0, 1);
        h_layer_energy_ratio[lIdx].Sumw2();
    }
}

inline std::shared_ptr<TH1D> buildHistoFromVector(
    const std::vector<double> &energyValues,
    const std::vector<double> &consgFactor)
{
    std::shared_ptr<TH1D> histo = std::make_shared<TH1D>("histo", "histoTitle", consgFactor.size(), energyValues[0], energyValues[energyValues.size() - 1]);
    for (auto idx = 1; idx <= histo->GetNbinsX(); ++idx)
        histo->SetBinContent(idx, consgFactor[idx - 1]);

    return histo;
}

inline void updateProcessStatus(const int evIdx, int &kStep, const int nevents)
{
    auto percentage = ((evIdx + 1) / (double)nevents) * 100;
    if (floor(percentage) != 0 && ((int)floor(percentage) % kStep) == 0)
    {
        std::cout << "\n"
                  << percentage << " %\t | \tProcessed " << evIdx + 1 << " events / " << nevents;
        kStep += 10;
    }
}

void buildAcceptance(
    const std::string accInputPath,
    const bool verbose,
    const std::vector<float> &logEBins,
    TFile &outFile,
    const std::string wd)
{

    //auto dmpch = aggregateEventsDmpChain(accInputPath,verbose);
    auto dmpch = aggregateEventsTChain(accInputPath, verbose);

    // Register Header container
    std::shared_ptr<DmpEvtHeader> evt_header = std::make_shared<DmpEvtHeader>();
    dmpch->SetBranchAddress("EventHeader", &evt_header);

    // Register SimuPrimaries container
    std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries = std::make_shared<DmpEvtSimuPrimaries>();
    dmpch->SetBranchAddress("DmpEvtSimuPrimaries", &simu_primaries);

    // Register BGO constainer
    std::shared_ptr<DmpEvtBgoHits> bgohits = std::make_shared<DmpEvtBgoHits>();
    dmpch->SetBranchAddress("DmpEvtBgoHits", &bgohits);

    // Register BGO REC constainer
    std::shared_ptr<DmpEvtBgoRec> bgorec = std::make_shared<DmpEvtBgoRec>();
    dmpch->SetBranchAddress("DmpEvtBgoRec", &bgorec);

    // Register STK container
    std::shared_ptr<TClonesArray> stkclusters = std::make_shared<TClonesArray>("DmpStkSiCluster");
    dmpch->SetBranchAddress("StkClusterCollection", &stkclusters);

    // Check if STK tracks collection exists
    bool fStkKalmanTracksFound = false;
    for (int brIdx = 0; brIdx < dmpch->GetListOfBranches()->GetEntries(); ++brIdx)
        if (strcmp(dmpch->GetListOfBranches()->At(brIdx)->GetName(), "StkKalmanTracks"))
        {
            fStkKalmanTracksFound = true;
            break;
        }

    // Register STK tracks collection
    std::shared_ptr<TClonesArray> stktracks = std::make_shared<TClonesArray>("DmpStkTrack", 200);
    if (fStkKalmanTracksFound)
        dmpch->SetBranchAddress("StkKalmanTracks", &stktracks);

    // Register PSD container
    std::shared_ptr<DmpEvtPsdHits> psdhits = std::make_shared<DmpEvtPsdHits>();
    dmpch->SetBranchAddress("DmpPsdHits", &psdhits);

    // Event loop
    auto nevents = dmpch->GetEntries();
    if (verbose)
        std::cout << "\n\nTotal number of events: " << nevents << "\n\n";

    // Acceptance - First-Cut histos
    TH1D h_geo_factor("h_geo_factor", "Energy Distribution of the geometric factor", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_incoming("h_incoming", "Energy Distribution of the incoming particles", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_trigger("h_trigger", "Energy Distribution of the triggered particles", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_gometric_cut("h_gometric_cut", "Energy Distribution - geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxElayer_cut("h_maxElayer_cut", "Energy Distribution - maxElayer cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxBarLayer_cut("h_maxBarLayer_cut", "Energy Distribution - maxBarLayer cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOTrackContainment_cut("h_BGOTrackContainment_cut", "Energy Distribution - BGOTrackContainment cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGO_fiducial_cut("h_BGO_fiducial_cut", "Energy Distibution - BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_nBarLayer13_cut("h_nBarLayer13_cut", "Energy Distribution - nBarLayer13 cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxRms_cut("h_maxRms_cut", "Energy Distribution - maxRms cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_track_selection_cut("h_track_selection_cut", "Energy Distribution - track selection cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_xtrl_cut("h_xtrl_cut", "Energy Distribution - xtrl cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_psd_charge_cut("h_psd_charge_cut", "Energy Distribution - psd charge cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_stk_charge_cut("h_stk_charge_cut", "Energy Distribution - stk charge cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_all_cut("h_all_cut", "Energy Distribution - All cut ", logEBins.size() - 1, &(logEBins[0]));

    TH1D h_geo_factor_w("h_geo_factor_w", "Energy Distribution of the geometric factor", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_incoming_w("h_incoming_w", "Energy Distribution of the incoming particles", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_trigger_w("h_trigger_w", "Energy Distribution of the triggered particles", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_gometric_cut_w("h_gometric_cut_w", "Energy Distribution - geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxElayer_cut_w("h_maxElayer_cut_w", "Energy Distribution - maxElayer cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxBarLayer_cut_w("h_maxBarLayer_cut_w", "Energy Distribution - maxBarLayer cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOTrackContainment_cut_w("h_BGOTrackContainment_cut_w", "Energy Distribution - BGOTrackContainment cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGO_fiducial_cut_w("h_BGO_fiducial_cut_w", "Energy Distibution - BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_nBarLayer13_cut_w("h_nBarLayer13_cut_w", "Energy Distribution - nBarLayer13 cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxRms_cut_w("h_maxRms_cut_w", "Energy Distribution - maxRms cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_track_selection_cut_w("h_track_selection_cut_w", "Energy Distribution - track selection cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_xtrl_cut_w("h_xtrl_cut_w", "Energy Distribution - xtrl cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_psd_charge_cut_w("h_psd_charge_cut_w", "Energy Distribution - psd charge cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_stk_charge_cut_w("h_stk_charge_cut_w", "Energy Distribution - stk charge cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_all_cut_w("h_all_cut_w", "Energy Distribution - All cut ", logEBins.size() - 1, &(logEBins[0]));

    // Acceptance - Cuts && Geometric Cut
    TH1D h_geometric_maxElayer_cut("h_geometric_maxElayer_cut", "Energy Distribution - maxElayer + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_maxBarLayer_cut("h_geometric_maxBarLayer_cut", "Energy Distribution - maxBarLayer + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_BGOTrackContainment_cut("h_geometric_BGOTrackContainment_cut", "Energy Distribution - BGOTrackContainment + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_BGO_fiducial_cut("h_geometric_BGO_fiducial_cut", "Energy Distibution - BGO fiducial + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_nBarLayer13_cut("h_geometric_nBarLayer13_cut", "Energy Distribution - nBarLayer13 + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_maxRms_cut("h_geometric_maxRms_cut", "Energy Distribution - maxRms + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_track_selection_cut("h_geometric_track_selection_cut", "Energy Distribution - track selection + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_xtrl_cut("h_geometric_xtrl_cut", "Energy Distribution - xtrl + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_psd_charge_cut("h_geometric_psd_charge_cut", "Energy Distribution - psd charge + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_stk_charge_cut("h_geometric_stk_charge_cut", "Energy Distribution - stk charge + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_all_cut("h_geometric_all_cut", "Energy Distribution - All + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));

    TH1D h_geometric_maxElayer_cut_w("h_geometric_maxElayer_cut_w", "Energy Distribution - maxElayer + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_maxBarLayer_cut_w("h_geometric_maxBarLayer_cut_w", "Energy Distribution - maxBarLayer + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_BGOTrackContainment_cut_w("h_geometric_BGOTrackContainment_cut_w", "Energy Distribution - BGOTrackContainment + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_BGO_fiducial_cut_w("h_geometric_BGO_fiducial_cut_w", "Energy Distibution - BGO fiducial + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_nBarLayer13_cut_w("h_geometric_nBarLayer13_cut_w", "Energy Distribution - nBarLayer13 + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_maxRms_cut_w("h_geometric_maxRms_cut_w", "Energy Distribution - maxRms + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_track_selection_cut_w("h_geometric_track_selection_cut_w", "Energy Distribution - track selection + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_xtrl_cut_w("h_geometric_xtrl_cut_w", "Energy Distribution - xtrl + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_psd_charge_cut_w("h_geometric_psd_charge_cut_w", "Energy Distribution - psd charge + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_stk_charge_cut_w("h_geometric_stk_charge_cut_w", "Energy Distribution - stk charge + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_all_cut_w("h_geometric_all_cut_w", "Energy Distribution - All + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));

    // Acceptance - Cuts && BGO fiducial volume cut
    TH1D h_BGOfiducial_nBarLayer13_cut("h_BGOfiducial_nBarLayer13_cut", "Energy Distribution - nBarLayer13 + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_maxRms_cut("h_BGOfiducial_maxRms_cut", "Energy Distribution - maxRms  + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_track_selection_cut("h_BGOfiducial_track_selection_cut", "Energy Distribution - track selection + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_xtrl_cut("h_BGOfiducial_xtrl_cut", "Energy Distribution - xtrl + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_psd_charge_cut("h_BGOfiducial_psd_charge_cut", "Energy Distribution - psd charge + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_stk_charge_cut("h_BGOfiducial_stk_charge_cut", "Energy Distribution - stk charge + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_all_cut("h_BGOfiducial_all_cut", "Energy Distribution - All + BGO fiducial cut ", logEBins.size() - 1, &(logEBins[0]));

    TH1D h_BGOfiducial_nBarLayer13_cut_w("h_BGOfiducial_nBarLayer13_cut_w", "Energy Distribution - nBarLayer13 + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_maxRms_cut_w("h_BGOfiducial_maxRms_cut_w", "Energy Distribution - maxRms  + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_track_selection_cut_w("h_BGOfiducial_track_selection_cut_w", "Energy Distribution - track selection + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_xtrl_cut_w("h_BGOfiducial_xtrl_cut_w", "Energy Distribution - xtrl + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_psd_charge_cut_w("h_BGOfiducial_psd_charge_cut_w", "Energy Distribution - psd charge + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_stk_charge_cut_w("h_BGOfiducial_stk_charge_cut_w", "Energy Distribution - stk charge + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_all_cut_w("h_BGOfiducial_all_cut_w", "Energy Distribution - All + BGO fiducial cut ", logEBins.size() - 1, &(logEBins[0]));

    // Analysis histos - simu and reco energy of incoming events
    TH1D h_BGOrec_E("h_BGOrec_E", "BGO Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOrec_E_corr("h_BGOrec_E_corr", "BGO Corrected Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_simu_energy("h_simu_energy", "Simu Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_energy_diff("h_energy_diff", "Simu vs Corrected Reco BGO energy", 50, -100, 100);

    // Analysis histos - simu and reco energy of triggered events
    TH1D h_triggered_BGOrec_E("h_triggered_BGOrec_E", "Triggered BGO Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_triggered_BGOrec_E_corr("h_triggered_BGOrec_E_corr", "Triggered BGO Corrected Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_triggered_simu_energy("h_triggered_simu_energy", "Triggered Simu Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_triggered_energy_diff("h_triggered_energy_diff", "Simu vs Corrected Triggered Reco BGO energy", 50, -100, 100);

    // Pre Geometric Cut
    // Top X and Y
    TH1D h_preGeo_BGOrec_topX_vs_realX("h_preGeo_BGOrec_topX_vs_realX", "Real X - BGOrec TOP X", 100, -100, 100);
    TH1D h_preGeo_BGOrec_topY_vs_realY("h_preGeo_BGOrec_topY_vs_realY", "Real Y - BGOrec TOP Y", 100, -100, 100);

    // Slope X and Y
    TH1D h_preGeo_real_slopeX("h_preGeo_real_slopeX", "Real Slope X", 1000, -90, 90);
    TH1D h_preGeo_real_slopeY("h_preGeo_real_slopeY", "Real Slope Y", 1000, -90, 90);
    TH1D h_preGeo_BGOrec_slopeX("h_preGeo_BGOrec_slopeX", "BGOrec Slope X", 1000, -90, 90);
    TH1D h_preGeo_BGOrec_slopeY("h_preGeo_BGOrec_slopeY", "BGOrec Slope Y", 1000, -90, 90);

    // Intercept X and Y
    TH1D h_preGeo_real_interceptX("h_preGeo_real_interceptX", "Real Intercept X", 500, -500, 500);
    TH1D h_preGeo_real_interceptY("h_preGeo_real_interceptY", "Real Intercept Y", 500, -500, 500);
    TH1D h_preGeo_BGOrec_interceptX("h_preGeo_BGOrec_interceptX", "BGOrec Intercept X", 500, -500, 500);
    TH1D h_preGeo_BGOrec_interceptY("h_preGeo_BGOrec_interceptY", "BGOrec Intercept Y", 500, -500, 500);

    // Top Maps
    TH2D h_preGeo_real_topMap("h_preGeo_real_topMap", "Real BGO TOP Map", 500, -500, 500, 500, -500, 500);
    TH2D h_preGeo_BGOreco_topMap("h_preGeo_BGOreco_topMap", "BGOreco TOP Map", 500, -500, 500, 500, -500, 500);

    // Bottom Maps
    TH2D h_preGeo_real_bottomMap("h_preGeo_real_bottomMap", "Real BGO BOTTOM Map", 500, -500, 500, 500, -500, 500);
    TH2D h_preGeo_BGOreco_bottomMap("h_preGeo_BGOreco_bottomMap", "BGOreco BOTTOM Map", 500, -500, 500, 500, -500, 500);

    // Map of events outside the "real" first BGO layer
    TH2D h_noBGOenergy_real_topMap("h_noBGOenergy_real_topMap", "Real BGO TOP Map", 500, -500, 500, 500, -500, 500);

    // After Geometric Cut
    // Top X and Y
    TH1D h_geo_BGOrec_topX_vs_realX("h_geo_BGOrec_topX_vs_realX", "Real X - BGOrec TOP X", 100, -100, 100);
    TH1D h_geo_BGOrec_topY_vs_realY("h_geo_BGOrec_topY_vs_realY", "Real Y - BGOrec TOP Y", 100, -100, 100);

    // Slope X and Y
    TH1D h_geo_real_slopeX("h_geo_real_slopeX", "Real Slope X", 1000, -90, 90);
    TH1D h_geo_real_slopeY("h_geo_real_slopeY", "Real Slope Y", 1000, -90, 90);
    TH1D h_geo_BGOrec_slopeX("h_geo_BGOrec_slopeX", "BGOrec Slope X", 1000, -90, 90);
    TH1D h_geo_BGOrec_slopeY("h_geo_BGOrec_slopeY", "BGOrec Slope Y", 1000, -90, 90);

    // Intercept X and Y
    TH1D h_geo_real_interceptX("h_geo_real_interceptX", "Real Intercept X", 500, -500, 500);
    TH1D h_geo_real_interceptY("h_geo_real_interceptY", "Real Intercept Y", 500, -500, 500);
    TH1D h_geo_BGOrec_interceptX("h_geo_BGOrec_interceptX", "BGOrec Intercept X", 500, -500, 500);
    TH1D h_geo_BGOrec_interceptY("h_geo_BGOrec_interceptY", "BGOrec Intercept Y", 500, -500, 500);

    // Top Maps
    TH2D h_geo_real_topMap("h_geo_real_topMap", "Real BGO TOP Map", 500, -500, 500, 500, -500, 500);
    TH2D h_geo_BGOreco_topMap("h_geo_BGOreco_topMap", "BGOreco TOP Map", 500, -500, 500, 500, -500, 500);

    // Bottom Maps
    TH2D h_geo_real_bottomMap("h_geo_real_bottomMap", "Real BGO BOTTOM Map", 500, -500, 500, 500, -500, 500);
    TH2D h_geo_BGOreco_bottomMap("h_geo_BGOreco_bottomMap", "BGOreco BOTTOM Map", 500, -500, 500, 500, -500, 500);

    // Ratio of layer energy respect to total BGO energy
    TH1D h_layer_max_energy_ratio("h_layer_max_energy_ratio", "Layer Energy Ratio", 100, 0, 1);
    std::vector<TH1D> h_layer_energy_ratio;
    init_BGO_histos(h_layer_energy_ratio);

    // XTRL histos
    auto xtrl_bins = LinearSpacedArray(0, 100, 1000);

    TH1D h_xtrl_energy_int("h_xtrl_energy_int", "Energy integrated XTRL distribution", xtrl_bins.size() - 1, &(xtrl_bins[0]));
    TH2D h_xtrl("h_xtrl", "XTRL energy Distribution", logEBins.size() - 1, &(logEBins[0]), xtrl_bins.size() - 1, &(xtrl_bins[0]));

    // STK charge histos
    TH1D h_chargeX("h_chargeX", "Charge distribution X", 1000, 0, 1000);
    TH1D h_chargeY("h_chargeY", "Charge distribution Y", 1000, 0, 1000);
    TH2D h_charge2D("h_charge2D", "STK charge", 1000, 0, 1000, 1000, 0, 1000);
    TH1D h_charge("h_charge", "Mean STK charge", 1000, 0, 1000);

    TH1D h_selected_chargeX("h_selected_chargeX", "Charge distribution X", 1000, 0, 1000);
    TH1D h_selected_chargeY("h_selected_chargeY", "Charge distribution Y", 1000, 0, 1000);
    TH2D h_selected_charge2D("h_selected_charge2D", "STK charge", 1000, 0, 1000, 1000, 0, 1000);
    TH1D h_selected_charge("h_selected_charge", "Mean STK charge", 1000, 0, 1000);

    // Proton background histos
    TH1D h_background_under_xtrl_cut("h_background_under_xtrl_cut", "Proton background - XTRL < cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_background_over_xtrl_cut("h_background_over_xtrl_cut", "Proton background - 20 < XTRL < 100", logEBins.size() - 1, &(logEBins[0]));

    // Sumw2 Acceptance - First-Cut histos
    h_geo_factor.Sumw2();
    h_incoming.Sumw2();
    h_trigger.Sumw2();
    h_gometric_cut.Sumw2();
    h_maxElayer_cut.Sumw2();
    h_maxBarLayer_cut.Sumw2();
    h_BGOTrackContainment_cut.Sumw2();
    h_BGO_fiducial_cut.Sumw2();
    h_nBarLayer13_cut.Sumw2();
    h_maxRms_cut.Sumw2();
    h_track_selection_cut.Sumw2();
    h_xtrl_cut.Sumw2();
    h_psd_charge_cut.Sumw2();
    h_stk_charge_cut.Sumw2();
    h_all_cut.Sumw2();

    h_geo_factor_w.Sumw2();
    h_incoming_w.Sumw2();
    h_trigger_w.Sumw2();
    h_gometric_cut_w.Sumw2();
    h_maxElayer_cut_w.Sumw2();
    h_maxBarLayer_cut_w.Sumw2();
    h_BGOTrackContainment_cut_w.Sumw2();
    h_BGO_fiducial_cut_w.Sumw2();
    h_nBarLayer13_cut_w.Sumw2();
    h_maxRms_cut_w.Sumw2();
    h_track_selection_cut_w.Sumw2();
    h_xtrl_cut_w.Sumw2();
    h_psd_charge_cut_w.Sumw2();
    h_stk_charge_cut_w.Sumw2();
    h_all_cut_w.Sumw2();

    // Sumw2 Acceptance - Cuts && Geometric Cut
    h_geometric_maxElayer_cut.Sumw2();
    h_geometric_maxBarLayer_cut.Sumw2();
    h_geometric_BGOTrackContainment_cut.Sumw2();
    h_geometric_BGO_fiducial_cut.Sumw2();
    h_geometric_nBarLayer13_cut.Sumw2();
    h_geometric_maxRms_cut.Sumw2();
    h_geometric_track_selection_cut.Sumw2();
    h_geometric_xtrl_cut.Sumw2();
    h_geometric_psd_charge_cut.Sumw2();
    h_geometric_stk_charge_cut.Sumw2();
    h_geometric_all_cut.Sumw2();

    h_geometric_maxElayer_cut_w.Sumw2();
    h_geometric_maxBarLayer_cut_w.Sumw2();
    h_geometric_BGOTrackContainment_cut_w.Sumw2();
    h_geometric_BGO_fiducial_cut_w.Sumw2();
    h_geometric_nBarLayer13_cut_w.Sumw2();
    h_geometric_maxRms_cut_w.Sumw2();
    h_geometric_track_selection_cut_w.Sumw2();
    h_geometric_xtrl_cut_w.Sumw2();
    h_geometric_psd_charge_cut_w.Sumw2();
    h_geometric_stk_charge_cut_w.Sumw2();
    h_geometric_all_cut_w.Sumw2();

    // Sumw2 Acceptance - Cuts && BGO fiducial volume cut
    h_BGOfiducial_nBarLayer13_cut.Sumw2();
    h_BGOfiducial_maxRms_cut.Sumw2();
    h_BGOfiducial_track_selection_cut.Sumw2();
    h_BGOfiducial_xtrl_cut.Sumw2();
    h_BGOfiducial_psd_charge_cut.Sumw2();
    h_BGOfiducial_stk_charge_cut.Sumw2();
    h_BGOfiducial_all_cut.Sumw2();

    h_BGOfiducial_nBarLayer13_cut_w.Sumw2();
    h_BGOfiducial_maxRms_cut_w.Sumw2();
    h_BGOfiducial_track_selection_cut_w.Sumw2();
    h_BGOfiducial_xtrl_cut_w.Sumw2();
    h_BGOfiducial_psd_charge_cut_w.Sumw2();
    h_BGOfiducial_stk_charge_cut_w.Sumw2();
    h_BGOfiducial_all_cut_w.Sumw2();

    // Sumw2 Analysis histos - simu and reco energy of incoming events
    h_BGOrec_E.Sumw2();
    h_BGOrec_E_corr.Sumw2();
    h_simu_energy.Sumw2();
    h_energy_diff.Sumw2();

    // Sumw2 Analysis histos - simu and reco energy of triggered events
    h_triggered_BGOrec_E.Sumw2();
    h_triggered_BGOrec_E_corr.Sumw2();
    h_triggered_simu_energy.Sumw2();
    h_triggered_energy_diff.Sumw2();

    // Sumw2 Analysis histos - preGeo
    h_preGeo_BGOrec_topX_vs_realX.Sumw2();
    h_preGeo_BGOrec_topY_vs_realY.Sumw2();
    h_preGeo_real_slopeX.Sumw2();
    h_preGeo_real_slopeY.Sumw2();
    h_preGeo_BGOrec_slopeX.Sumw2();
    h_preGeo_BGOrec_slopeY.Sumw2();
    h_preGeo_real_interceptX.Sumw2();
    h_preGeo_real_interceptY.Sumw2();
    h_preGeo_BGOrec_interceptX.Sumw2();
    h_preGeo_BGOrec_interceptY.Sumw2();
    h_preGeo_real_topMap.Sumw2();
    h_preGeo_BGOreco_topMap.Sumw2();
    h_preGeo_real_bottomMap.Sumw2();
    h_preGeo_BGOreco_bottomMap.Sumw2();
    h_noBGOenergy_real_topMap.Sumw2();

    // Sumw2 Analysis histos - Geo
    h_geo_BGOrec_topX_vs_realX.Sumw2();
    h_geo_BGOrec_topY_vs_realY.Sumw2();
    h_geo_real_slopeX.Sumw2();
    h_geo_real_slopeY.Sumw2();
    h_geo_BGOrec_slopeX.Sumw2();
    h_geo_BGOrec_slopeY.Sumw2();
    h_geo_real_interceptX.Sumw2();
    h_geo_real_interceptY.Sumw2();
    h_geo_BGOrec_interceptX.Sumw2();
    h_geo_BGOrec_interceptY.Sumw2();
    h_geo_real_topMap.Sumw2();
    h_geo_BGOreco_topMap.Sumw2();
    h_geo_real_bottomMap.Sumw2();
    h_geo_BGOreco_bottomMap.Sumw2();
    h_layer_max_energy_ratio.Sumw2();

    // Sumw2 XTRL histos
    h_xtrl_energy_int.Sumw2();
    h_xtrl.Sumw2();

    // Sumw2 STK charge histos
    h_chargeX.Sumw2();
    h_chargeY.Sumw2();
    h_charge.Sumw2();
    h_charge2D.Sumw2();

    h_selected_chargeX.Sumw2();
    h_selected_chargeY.Sumw2();
    h_selected_charge.Sumw2();
    h_selected_charge2D.Sumw2();
    
    // Proton background histos
    h_background_under_xtrl_cut.Sumw2();
    h_background_over_xtrl_cut.Sumw2();

    // Create and load acceptance events cuts from config file
    
    // Create acceptance cuts struct
    cuts_conf acceptance_cuts;
    // Create active cuts struct
    data_active_cuts active_cuts;
    // Create ancillary cuts struct
    mc_ancillary_cuts ancillary_cuts;

    // Load structs reading config file
    load_acceptance_struct(
        acceptance_cuts, 
        active_cuts, 
        ancillary_cuts, 
        wd);

    // Read dataSets connfig file
    data_set_conf input_sets;
    load_input_dsets_config(input_sets, wd);

    double _GeV = 0.001;
    int kStep = 10;

    for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx)
    {
        // Get chain event
        dmpch->GetEvent(evIdx);

        // Event printout
        if (verbose)
            updateProcessStatus(evIdx, kStep, nevents);

        // Get event total energy
        double bgoTotalE_raw = bgorec->GetTotalEnergy(); // Energy in MeV - not corrected
        double bgoTotalE = bgorec->GetElectronEcor();    // Returns corrected energy assuming this was an electron (MeV)
        double simuEnergy = simu_primaries->pvpart_ekin; //Energy of simu primaries particle (MeV)

        double energy_w = pow(simuEnergy * _GeV, -1.7);        // Set the energy histo weight -> Used to reweight 1/E energy histos

        // Fill the energy histos
        h_BGOrec_E.Fill(bgoTotalE_raw * _GeV);
        h_BGOrec_E_corr.Fill(bgoTotalE * _GeV);
        h_simu_energy.Fill(simuEnergy * _GeV);
        h_energy_diff.Fill((simuEnergy - bgoTotalE) * _GeV);

        // Don't accept events outside the selected energy window
        if (simuEnergy * _GeV < acceptance_cuts.min_event_energy || simuEnergy * _GeV > acceptance_cuts.max_event_energy)
            continue;

        if (geometric_cut(simu_primaries))
        {
            h_geo_factor.Fill(simuEnergy * _GeV);
            h_geo_factor_w.Fill(simuEnergy * _GeV, energy_w);
        }

        // Read trigger status
        // For MC events triggers 1 and 2 are always disabled
        bool unbiased_tr = evt_header->GeneratedTrigger(0) && evt_header->EnabledTrigger(0);
        bool mip1_tr = evt_header->GeneratedTrigger(1);
        bool mip2_tr = evt_header->GeneratedTrigger(2);
        bool HET_tr = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
        bool LET_tr = evt_header->GeneratedTrigger(4) && evt_header->EnabledTrigger(4);
        bool MIP_tr = mip1_tr || mip2_tr;

        bool general_trigger = MIP_tr || HET_tr || LET_tr;

        // Check if the event has been triggered or not
        if (general_trigger)
        {
            if (checkBGOreco(bgorec, simu_primaries))
            {
                h_trigger.Fill(simuEnergy * _GeV);
                h_trigger_w.Fill(simuEnergy * _GeV, energy_w);
                h_incoming.Fill(simuEnergy * _GeV);
                h_incoming_w.Fill(simuEnergy * _GeV, energy_w);
            }
            else
                continue;
        }
        else
        {
            // Increase the incoming events
            h_incoming.Fill(simuEnergy * _GeV);
            h_incoming.Fill(simuEnergy * _GeV, energy_w);

            // Evaluate the position on the First BGO layer for the non - triggered events
            evaluateTopBottomPosition(
                simu_primaries,
                bgorec,
                h_preGeo_BGOrec_topX_vs_realX,
                h_preGeo_BGOrec_topY_vs_realY,
                h_preGeo_real_slopeX,
                h_preGeo_real_slopeY,
                h_preGeo_BGOrec_slopeX,
                h_preGeo_BGOrec_slopeY,
                h_preGeo_real_interceptX,
                h_preGeo_real_interceptY,
                h_preGeo_BGOrec_interceptX,
                h_preGeo_BGOrec_interceptY,
                h_preGeo_real_topMap,
                h_preGeo_BGOreco_topMap,
                h_preGeo_real_bottomMap,
                h_preGeo_BGOreco_bottomMap);

            continue;
        }

        best_track event_best_track;

        // Fill energy histos
        h_triggered_BGOrec_E.Fill(bgoTotalE_raw * _GeV);
        h_triggered_BGOrec_E_corr.Fill(bgoTotalE * _GeV);
        h_triggered_simu_energy.Fill(simuEnergy * _GeV);
        h_triggered_energy_diff.Fill((simuEnergy - bgoTotalE) * _GeV);

        // Load BGO event class
        DmpBgoContainer bgoVault(DAMPE_bgo_nLayers);
        bgoVault.scanBGOHits(
            bgohits,
            bgoTotalE,
            DAMPE_bgo_nLayers);

        // evaluate the energy raio on each single layer of the BGO
        evaluateEnergyRatio(
            bgorec,
            acceptance_cuts,
            bgoTotalE,
            h_layer_max_energy_ratio,
            h_layer_energy_ratio);

        auto filter_geometric_cut = false;
        auto filter_BGO_fiducial_cut = false;
        auto filter_BGO_fiducial_maxElayer_cut = false;
        auto filter_BGO_fiducial_maxBarLayer_cut = false;
        auto filter_BGO_fiducial_BGOTrackContainment_cut = false;
        auto filter_nBarLayer13_cut = false;
        auto filter_maxRms_cut = false;
        auto filter_track_selection_cut = false;
        auto filter_xtrl_cut = false;
        auto filter_psd_charge_cut = false;
        auto filter_stk_charge_cut = false;
        auto filter_all_cut = true;

        // Cut check...

        // **** Geometric cut ****
        filter_geometric_cut = geometric_cut(simu_primaries);

        // **** BGO Fiducial Volume ****
        if (active_cuts.BGO_fiducial)
        {
            // maxElayer_cut
            filter_BGO_fiducial_maxElayer_cut = maxElayer_cut(
                bgorec,
                acceptance_cuts,
                bgoTotalE);

            // maxBarLayer_cut
            filter_BGO_fiducial_maxBarLayer_cut = maxBarLayer_cut(
                bgoVault.GetLayerBarNumber(),
                bgoVault.GetiMaxLayer(),
                bgoVault.GetIdxBarMaxLayer());

            // BGOTrackContainment_cut
            filter_BGO_fiducial_BGOTrackContainment_cut = BGOTrackContainment_cut(
                bgorec,
                acceptance_cuts);

            filter_BGO_fiducial_cut = filter_BGO_fiducial_maxElayer_cut && filter_BGO_fiducial_maxBarLayer_cut && filter_BGO_fiducial_BGOTrackContainment_cut;
            filter_all_cut *= filter_BGO_fiducial_cut;
        }

        // **** nBarLayer13 cut ****
        if (active_cuts.nBarLayer13)
        {
            filter_nBarLayer13_cut = nBarLayer13_cut(
                bgohits,
                bgoVault.GetSingleLayerBarNumber(13),
                bgoTotalE);
            filter_all_cut *= filter_nBarLayer13_cut;
        }

        // **** maxRms cut ****
        if (active_cuts.maxRms)
        {
            filter_maxRms_cut = maxRms_cut(
                bgoVault.GetLayerBarNumber(),
                bgoVault.GetRmsLayer(),
                bgoTotalE,
                acceptance_cuts);
            filter_all_cut *= filter_maxRms_cut;
        }

        // **** track_selection cut ****
        if (active_cuts.track_selection)
        {
            filter_track_selection_cut =
                track_selection_cut(
                    bgorec,
                    bgohits,
                    stkclusters,
                    stktracks,
                    acceptance_cuts,
                    event_best_track);
            filter_all_cut *= filter_track_selection_cut;
        }

        // **** xtrl cut ****
        if (active_cuts.xtrl)
        {
            filter_xtrl_cut = xtrl_cut(
                bgoVault.GetSumRMS(),
                bgoVault.GetFracLayer(),
                acceptance_cuts,
                simuEnergy,
                h_xtrl_energy_int,
                h_xtrl);
            filter_all_cut *= filter_xtrl_cut;
        }

        // **** psd_charge cut ****
        if (active_cuts.psd_charge)
        {
            if (active_cuts.track_selection)
            {
                if (filter_track_selection_cut)
                {
                    filter_psd_charge_cut = psd_charge_cut(
                        psdhits,
                        bgorec,
                        acceptance_cuts,
                        event_best_track);
                    filter_all_cut *= filter_psd_charge_cut;
                }
            }
        }

        // **** stk_charge cut ****
        if (active_cuts.stk_charge)
        {   
            if (active_cuts.track_selection)
            {
                if (filter_track_selection_cut)
                {
                    // Fill charge histos
                    fillChargeHistos(
                        h_chargeX, 
                        h_chargeY,
                        h_charge,
                        h_charge2D,
                        event_best_track,
                        stkclusters);
            
                    // Charge cut
                    filter_stk_charge_cut = stk_charge_cut(
                        event_best_track,
                        stkclusters,
                        acceptance_cuts);
                    filter_all_cut *= filter_stk_charge_cut;
                }
            }
        }

        // **** ANCILLARY CUTS ****

        // **** compute proton background ****
        if (ancillary_cuts.compute_proton_background)
            compute_proton_background(
                bgoVault.GetSumRMS(),
                bgoVault.GetFracLayer(),
                acceptance_cuts,
                simuEnergy,
                h_background_under_xtrl_cut,
                h_background_over_xtrl_cut);


        // Fill cuts histos

        // Fill geometric cut histos
        if (filter_geometric_cut)
        {
            h_gometric_cut.Fill(simuEnergy * _GeV);
            h_gometric_cut_w.Fill(simuEnergy * _GeV, energy_w);

            // Evaluate the position on the First BGO layer (after geometric cut)
            evaluateTopBottomPosition(
                simu_primaries,
                bgorec,
                h_geo_BGOrec_topX_vs_realX,
                h_geo_BGOrec_topY_vs_realY,
                h_geo_real_slopeX,
                h_geo_real_slopeY,
                h_geo_BGOrec_slopeX,
                h_geo_BGOrec_slopeY,
                h_geo_real_interceptX,
                h_geo_real_interceptY,
                h_geo_BGOrec_interceptX,
                h_geo_BGOrec_interceptY,
                h_geo_real_topMap,
                h_geo_BGOreco_topMap,
                h_geo_real_bottomMap,
                h_geo_BGOreco_bottomMap);

            // Geometric cut && maxElayer cut
            if (filter_BGO_fiducial_maxElayer_cut)
            {
                h_geometric_maxElayer_cut.Fill(simuEnergy * _GeV);
                h_geometric_maxElayer_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // Geometric cut && maxBarLayer cut
            if (filter_BGO_fiducial_maxBarLayer_cut)
            {
                h_geometric_maxBarLayer_cut.Fill(simuEnergy * _GeV);
                h_geometric_maxBarLayer_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // Geometric cut && BGOTrackContainment cut
            if (filter_BGO_fiducial_BGOTrackContainment_cut)
            {
                h_geometric_BGOTrackContainment_cut.Fill(simuEnergy * _GeV);
                h_geometric_BGOTrackContainment_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // Geometric cut && BGO fiducial cut
            if (filter_BGO_fiducial_cut)
            {
                h_geometric_BGO_fiducial_cut.Fill(simuEnergy * _GeV);
                h_geometric_BGO_fiducial_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // Geometric cut && nBarLayer13 cut
            if (filter_nBarLayer13_cut)
            {
                h_geometric_nBarLayer13_cut.Fill(simuEnergy * _GeV);
                h_geometric_nBarLayer13_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // Geometric cut && maxRms cut
            if (filter_maxRms_cut)
            {
                h_geometric_maxRms_cut.Fill(simuEnergy * _GeV);
                h_geometric_maxRms_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // Geometric cut && track selection cut
            if (filter_track_selection_cut)
            {
                h_geometric_track_selection_cut.Fill(simuEnergy * _GeV);
                h_geometric_track_selection_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // Geometric cut && XTRL cut
            if (filter_xtrl_cut)
            {
                h_geometric_xtrl_cut.Fill(simuEnergy * _GeV);
                h_geometric_xtrl_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // Geometric cut && PSD charge cut
            if (filter_psd_charge_cut)
            {
                h_geometric_psd_charge_cut.Fill(simuEnergy * _GeV);
                h_geometric_psd_charge_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // Geometric cut && STK charge cut
            if (filter_stk_charge_cut)
            {
                h_geometric_stk_charge_cut.Fill(simuEnergy * _GeV);
                h_geometric_stk_charge_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // Geometric cut and all cuts
            if (filter_all_cut)
            {
                h_geometric_all_cut.Fill(simuEnergy * _GeV);
                h_geometric_all_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }
        }

        // Fill BGO_fiducial_maxElayer cut histo
        if (filter_BGO_fiducial_maxElayer_cut)
        {
            h_maxElayer_cut.Fill(simuEnergy * _GeV);
            h_maxElayer_cut_w.Fill(simuEnergy * _GeV, energy_w);
        }

        // Fill BGO_fiducial_maxBarLayer cut histo
        if (filter_BGO_fiducial_maxBarLayer_cut)
        {
            h_maxBarLayer_cut.Fill(simuEnergy * _GeV);
            h_maxBarLayer_cut_w.Fill(simuEnergy * _GeV, energy_w);
        }

        // Fill BGO_fiducial_BGOTrackContainment cut histo
        if (filter_BGO_fiducial_BGOTrackContainment_cut)
        {
            h_BGOTrackContainment_cut.Fill(simuEnergy * _GeV);
            h_BGOTrackContainment_cut_w.Fill(simuEnergy * _GeV, energy_w);
        }

        // Fill BGO fiducial volume cut
        if (filter_BGO_fiducial_cut)
        {
            h_BGO_fiducial_cut.Fill(simuEnergy * _GeV);
            h_BGO_fiducial_cut_w.Fill(simuEnergy * _GeV, energy_w);

            // BGO fiducial cut && nBarLayer13 cut
            if (filter_nBarLayer13_cut)
            {
                h_BGOfiducial_nBarLayer13_cut.Fill(simuEnergy * _GeV);
                h_BGOfiducial_nBarLayer13_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // BGO fiducial cut && maxRms cut
            if (filter_maxRms_cut)
            {
                h_BGOfiducial_maxRms_cut.Fill(simuEnergy * _GeV);
                h_BGOfiducial_maxRms_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // BGO fiducial cut && track selection cut
            if (filter_track_selection_cut)
            {
                h_BGOfiducial_track_selection_cut.Fill(simuEnergy * _GeV);
                h_BGOfiducial_track_selection_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // BGO fiducial cut && XTRL cut
            if (filter_xtrl_cut)
            {
                h_BGOfiducial_xtrl_cut.Fill(simuEnergy * _GeV);
                h_BGOfiducial_xtrl_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // BGO fiducial cut && PSD charge cut
            if (filter_psd_charge_cut)
            {
                h_BGOfiducial_psd_charge_cut.Fill(simuEnergy * _GeV);
                h_BGOfiducial_psd_charge_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // BGO fiducial cut && STK charge cut
            if (filter_stk_charge_cut)
            {
                h_BGOfiducial_stk_charge_cut.Fill(simuEnergy * _GeV);
                h_BGOfiducial_stk_charge_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }

            // BGO fiducial cut && all cut
            if (filter_all_cut)
            {
                h_BGOfiducial_all_cut.Fill(simuEnergy * _GeV);
                h_BGOfiducial_all_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }
        }

        // Fill nBarLayer13 cut histo
        if (filter_nBarLayer13_cut)
        {
            h_nBarLayer13_cut.Fill(simuEnergy * _GeV);
            h_nBarLayer13_cut_w.Fill(simuEnergy * _GeV, energy_w);
        }

        // Fill maxRms cut histo
        if (filter_maxRms_cut)
        {
            h_maxRms_cut.Fill(simuEnergy * _GeV);
            h_maxRms_cut_w.Fill(simuEnergy * _GeV, energy_w);
        }

        // Fill track selection cut histo
        if (filter_track_selection_cut)
        {
            h_track_selection_cut.Fill(simuEnergy * _GeV);
            h_track_selection_cut_w.Fill(simuEnergy * _GeV, energy_w);
        }

        // Fill XTRL cut histo
        if (filter_xtrl_cut)
        {
            h_xtrl_cut.Fill(simuEnergy * _GeV);
            h_xtrl_cut_w.Fill(simuEnergy * _GeV, energy_w);
        }

        // Fill PSD charge cut histo
        if (filter_psd_charge_cut)
        {
            h_psd_charge_cut.Fill(simuEnergy * _GeV);
            h_psd_charge_cut_w.Fill(simuEnergy * _GeV, energy_w);
        }

        // Fill STK charge histo
        if (filter_stk_charge_cut)
        {   
            if (filter_track_selection_cut)
            {
                // Fill selected STK charge histos
                fillChargeHistos(
                    h_selected_chargeX, 
                    h_selected_chargeY,
                    h_selected_charge,
                    h_selected_charge2D,
                    event_best_track,
                    stkclusters);

                h_stk_charge_cut.Fill(simuEnergy * _GeV);
                h_stk_charge_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }
        }

        // Fill all cut histo
        if (active_cuts.nActiveCuts)
            if (filter_all_cut)
            {
                h_all_cut.Fill(simuEnergy * _GeV);
                h_all_cut_w.Fill(simuEnergy * _GeV, energy_w);
            }
    }

    if (verbose)
    {
        std::cout << "\n\n ****** \n\n";
        std::cout << "generated events in good energy range: " << h_incoming.GetEntries() << std::endl;
        std::cout << "triggered events: " << h_trigger.GetEntries() << std::endl;

        auto refEntries = h_trigger.GetEntries();

        if (h_geo_factor.GetEntries())
            std::cout << "geometric filtered events (before trigger): " << h_geo_factor.GetEntries() << "/" << h_incoming.GetEntries() << " | statistic efficiency: " << static_cast<double>(h_geo_factor.GetEntries()) / h_incoming.GetEntries() << std::endl;

        if (h_gometric_cut.GetEntries())
            std::cout << "geometric filtered events (after trigger): " << h_gometric_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_gometric_cut.GetEntries()) / refEntries << std::endl;

        if (h_BGO_fiducial_cut.GetEntries())
            std::cout << "BGO fiducial filtered events: " << h_BGO_fiducial_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_BGO_fiducial_cut.GetEntries()) / refEntries << std::endl;

        if (h_nBarLayer13_cut.GetEntries())
            std::cout << "nBarLayer13 filtered events: " << h_nBarLayer13_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_nBarLayer13_cut.GetEntries()) / refEntries << std::endl;

        if (h_maxRms_cut.GetEntries())
            std::cout << "maxRms filtered events: " << h_maxRms_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_maxRms_cut.GetEntries()) / refEntries << std::endl;

        if (h_track_selection_cut.GetEntries())
            std::cout << "track_selection filtered events: " << h_track_selection_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_track_selection_cut.GetEntries()) / refEntries << std::endl;

        if (h_xtrl_cut.GetEntries())
            std::cout << "xtrl filtered events: " << h_xtrl_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_xtrl_cut.GetEntries()) / refEntries << std::endl;

        if (h_psd_charge_cut.GetEntries())
            std::cout << "psd_charge filtered events: " << h_psd_charge_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_psd_charge_cut.GetEntries()) / refEntries << std::endl;

        if (h_stk_charge_cut.GetEntries())
            std::cout << "stk_charge filtered events: " << h_stk_charge_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_stk_charge_cut.GetEntries()) / refEntries << std::endl;

        if (h_all_cut.GetEntries())
            std::cout << "all-cuts filtered events: " << h_all_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_all_cut.GetEntries()) / refEntries;

        std::cout << "\n\n ****** \n\n";
    }

    double genSurface = 4 * TMath::Pi() * pow(acceptance_cuts.vertex_radius, 2) / 2;
    double scaleFactor = TMath::Pi() * genSurface;

    // Building acceptance histos
    auto h_acceptance_geometric_factor = static_cast<TH1D *>(h_geo_factor.Clone("h_acceptance_geometric_factor"));
    auto h_acceptance_gometric_cut = static_cast<TH1D *>(h_gometric_cut.Clone("h_acceptance_gometric_cut"));
    auto h_acceptance_maxElayer_cut = static_cast<TH1D *>(h_maxElayer_cut.Clone("h_acceptance_maxElayer_cut"));
    auto h_acceptance_maxBarLayer_cut = static_cast<TH1D *>(h_maxBarLayer_cut.Clone("h_acceptance_maxBarLayer_cut"));
    auto h_acceptance_BGOTrackContainment_cut = static_cast<TH1D *>(h_BGOTrackContainment_cut.Clone("h_acceptance_BGOTrackContainment_cut"));
    auto h_acceptance_BGO_fiducial_cut = static_cast<TH1D *>(h_BGO_fiducial_cut.Clone("h_acceptance_BGO_fiducial_cut"));
    auto h_acceptance_BGO_fiducial_nBarLayer13_cut = static_cast<TH1D *>(h_BGOfiducial_nBarLayer13_cut.Clone("h_acceptance_BGO_fiducial_nBarLayer13_cut"));
    auto h_acceptance_BGO_fiducial_maxRms_cut = static_cast<TH1D *>(h_BGOfiducial_maxRms_cut.Clone("h_acceptance_BGO_fiducial_maxRms_cut"));
    auto h_acceptance_BGO_fiducial_track_selection_cut = static_cast<TH1D *>(h_BGOfiducial_track_selection_cut.Clone("h_acceptance_BGO_fiducial_track_selection_cut"));
    auto h_acceptance_BGO_fiducial_xtrl_cut = static_cast<TH1D *>(h_BGOfiducial_xtrl_cut.Clone("h_acceptance_BGO_fiducial_xtrl_cut"));
    auto h_acceptance_BGO_fiducial_psd_charge_cut = static_cast<TH1D *>(h_BGOfiducial_psd_charge_cut.Clone("h_acceptance_BGO_fiducial_psd_charge_cut"));
    auto h_acceptance_nBarLayer13_cut = static_cast<TH1D *>(h_nBarLayer13_cut.Clone("h_acceptance_nBarLayer13_cut"));
    auto h_acceptance_maxRms_cut = static_cast<TH1D *>(h_maxRms_cut.Clone("h_acceptance_maxRms_cut"));
    auto h_acceptance_track_selection_cut = static_cast<TH1D *>(h_track_selection_cut.Clone("h_acceptance_track_selection_cut"));
    auto h_acceptance_xtrl_cut = static_cast<TH1D *>(h_xtrl_cut.Clone("h_acceptance_xtrl_cut"));
    auto h_acceptance_psd_charge_cut = static_cast<TH1D *>(h_psd_charge_cut.Clone("h_acceptance_psd_charge_cut"));
    auto h_acceptance_stk_charge_cut = static_cast<TH1D *>(h_stk_charge_cut.Clone("h_acceptance_stk_charge_cut"));
    auto h_acceptance_all_cut = static_cast<TH1D *>(h_all_cut.Clone("h_acceptance_all_cut"));

    h_acceptance_geometric_factor->Divide(&h_incoming);
    h_acceptance_gometric_cut->Divide(&h_incoming);
    h_acceptance_maxElayer_cut->Divide(&h_incoming);
    h_acceptance_maxBarLayer_cut->Divide(&h_incoming);
    h_acceptance_BGOTrackContainment_cut->Divide(&h_incoming);
    h_acceptance_BGO_fiducial_cut->Divide(&h_incoming);
    h_acceptance_BGO_fiducial_nBarLayer13_cut->Divide(&h_incoming);
    h_acceptance_BGO_fiducial_maxRms_cut->Divide(&h_incoming);
    h_acceptance_BGO_fiducial_track_selection_cut->Divide(&h_incoming);
    h_acceptance_BGO_fiducial_xtrl_cut->Divide(&h_incoming);
    h_acceptance_BGO_fiducial_psd_charge_cut->Divide(&h_incoming);
    h_acceptance_nBarLayer13_cut->Divide(&h_incoming);
    h_acceptance_maxRms_cut->Divide(&h_incoming);
    h_acceptance_track_selection_cut->Divide(&h_incoming);
    h_acceptance_xtrl_cut->Divide(&h_incoming);
    h_acceptance_psd_charge_cut->Divide(&h_incoming);
    h_acceptance_stk_charge_cut->Divide(&h_incoming);
    h_acceptance_all_cut->Divide(&h_incoming);

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
    h_acceptance_stk_charge_cut->Scale(scaleFactor);
    h_acceptance_all_cut->Scale(scaleFactor);

    // Builing vectors
    std::vector<double> energyValues(h_incoming.GetXaxis()->GetNbins(), 0);

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
    std::vector<double> acceptanceValues_stk_charge_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_all_cut(energyValues.size(), 0);

    //Building histo errors on energy and
    std::vector<double> acceptanceError_geometric_factor(h_incoming.GetXaxis()->GetNbins(), 0);
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
    std::vector<double> acceptanceError_stk_charge_cut(acceptanceError_geometric_factor.size(), 0);
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
        acceptanceValues_stk_charge_cut[index] = h_acceptance_stk_charge_cut->GetBinContent(index + 1);
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
        acceptanceError_stk_charge_cut[index] = h_acceptance_stk_charge_cut->GetBinError(index + 1);
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
    TGraphErrors gr_acceptance_stk_charge_cut(energyValues.size(), &energyValues[0], &acceptanceValues_stk_charge_cut[0], &(energyError[0]), &(acceptanceError_stk_charge_cut[0]));
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
    gr_acceptance_stk_charge_cut.SetName("gr_acceptance_stk_charge_cut");
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
    gr_acceptance_stk_charge_cut.SetTitle("Acceptance - STK charge selection cut");
    gr_acceptance_all_cut.SetTitle("Acceptance - all cut");

    // Write histos to file
    // Acceptance - First-Cut histos
    h_geo_factor.Write();
    h_incoming.Write();
    h_trigger.Write();
    h_gometric_cut.Write();
    h_maxElayer_cut.Write();
    h_maxBarLayer_cut.Write();
    h_BGOTrackContainment_cut.Write();
    h_BGO_fiducial_cut.Write();
    h_nBarLayer13_cut.Write();
    h_maxRms_cut.Write();
    h_track_selection_cut.Write();
    h_xtrl_cut.Write();
    h_psd_charge_cut.Write();
    h_stk_charge_cut.Write();
    h_all_cut.Write();

    // Acceptance - Cuts && Geometric Cut
    h_geometric_maxElayer_cut.Write();
    h_geometric_maxBarLayer_cut.Write();
    h_geometric_BGOTrackContainment_cut.Write();
    h_geometric_BGO_fiducial_cut.Write();
    h_geometric_nBarLayer13_cut.Write();
    h_geometric_maxRms_cut.Write();
    h_geometric_track_selection_cut.Write();
    h_geometric_xtrl_cut.Write();
    h_geometric_psd_charge_cut.Write();
    h_geometric_stk_charge_cut.Write();
    h_geometric_all_cut.Write();
    // Acceptance - Cuts && BGO fiducial volume cut
    h_BGOfiducial_nBarLayer13_cut.Write();
    h_BGOfiducial_maxRms_cut.Write();
    h_BGOfiducial_track_selection_cut.Write();
    h_BGOfiducial_xtrl_cut.Write();
    h_BGOfiducial_psd_charge_cut.Write();
    h_BGOfiducial_stk_charge_cut.Write();
    h_BGOfiducial_all_cut.Write();

    // Writing reweighted histos
    auto reweightedHistoDir = outFile.mkdir("reweightedHistoDir");
    reweightedHistoDir->cd();

    // Acceptance - First-Cut histos
    h_geo_factor_w.Write();
    h_incoming_w.Write();
    h_trigger_w.Write();
    h_gometric_cut_w.Write();
    h_maxElayer_cut_w.Write();
    h_maxBarLayer_cut_w.Write();
    h_BGOTrackContainment_cut_w.Write();
    h_BGO_fiducial_cut_w.Write();
    h_nBarLayer13_cut_w.Write();
    h_maxRms_cut_w.Write();
    h_track_selection_cut_w.Write();
    h_xtrl_cut_w.Write();
    h_psd_charge_cut_w.Write();
    h_stk_charge_cut_w.Write();
    h_all_cut_w.Write();

    // Acceptance - Cuts && Geometric Cut
    h_geometric_maxElayer_cut_w.Write();
    h_geometric_maxBarLayer_cut_w.Write();
    h_geometric_BGOTrackContainment_cut_w.Write();
    h_geometric_BGO_fiducial_cut_w.Write();
    h_geometric_nBarLayer13_cut_w.Write();
    h_geometric_maxRms_cut_w.Write();
    h_geometric_track_selection_cut_w.Write();
    h_geometric_xtrl_cut_w.Write();
    h_geometric_psd_charge_cut_w.Write();
    h_geometric_stk_charge_cut_w.Write();
    h_geometric_all_cut_w.Write();
    // Acceptance - Cuts && BGO fiducial volume cut
    h_BGOfiducial_nBarLayer13_cut_w.Write();
    h_BGOfiducial_maxRms_cut_w.Write();
    h_BGOfiducial_track_selection_cut_w.Write();
    h_BGOfiducial_xtrl_cut_w.Write();
    h_BGOfiducial_psd_charge_cut_w.Write();
    h_BGOfiducial_stk_charge_cut_w.Write();
    h_BGOfiducial_all_cut_w.Write();

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
    h_acceptance_stk_charge_cut->Write();
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
    gr_acceptance_stk_charge_cut.Write();
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
    std::shared_ptr<TEfficiency> trigger_efficiency;
    std::shared_ptr<TEfficiency> tr_eff_gometric_cut;
    std::shared_ptr<TEfficiency> tr_eff_maxElayer_cut;
    std::shared_ptr<TEfficiency> tr_eff_maxBarLayer_cut;
    std::shared_ptr<TEfficiency> tr_eff_BGOTrackContainment_cut;
    std::shared_ptr<TEfficiency> tr_eff_BGO_fiducial_cut;
    std::shared_ptr<TEfficiency> tr_eff_nBarLayer13_cut;
    std::shared_ptr<TEfficiency> tr_eff_maxRms_cut;
    std::shared_ptr<TEfficiency> tr_eff_track_selection_cut;
    std::shared_ptr<TEfficiency> tr_eff_xtrl_cut;
    std::shared_ptr<TEfficiency> tr_eff_psd_charge_cut;
    std::shared_ptr<TEfficiency> tr_eff_stk_charge_cut;
    std::shared_ptr<TEfficiency> tr_eff_all_cut;

    if (TEfficiency::CheckConsistency(h_gometric_cut, h_gometric_cut))
        trigger_efficiency = std::make_shared<TEfficiency>(h_gometric_cut, h_gometric_cut);

    if (TEfficiency::CheckConsistency(h_gometric_cut, h_trigger))
        tr_eff_gometric_cut = std::make_shared<TEfficiency>(h_gometric_cut, h_trigger);

    if (TEfficiency::CheckConsistency(h_maxElayer_cut, h_trigger))
        tr_eff_maxElayer_cut = std::make_shared<TEfficiency>(h_maxElayer_cut, h_trigger);

    if (TEfficiency::CheckConsistency(h_maxBarLayer_cut, h_trigger))
        tr_eff_maxBarLayer_cut = std::make_shared<TEfficiency>(h_maxBarLayer_cut, h_trigger);

    if (TEfficiency::CheckConsistency(h_BGOTrackContainment_cut, h_trigger))
        tr_eff_BGOTrackContainment_cut = std::make_shared<TEfficiency>(h_BGOTrackContainment_cut, h_trigger);

    if (TEfficiency::CheckConsistency(h_BGO_fiducial_cut, h_trigger))
        tr_eff_BGO_fiducial_cut = std::make_shared<TEfficiency>(h_BGO_fiducial_cut, h_trigger);

    if (TEfficiency::CheckConsistency(h_nBarLayer13_cut, h_trigger))
        tr_eff_nBarLayer13_cut = std::make_shared<TEfficiency>(h_nBarLayer13_cut, h_trigger);

    if (TEfficiency::CheckConsistency(h_maxRms_cut, h_trigger))
        tr_eff_maxRms_cut = std::make_shared<TEfficiency>(h_maxRms_cut, h_trigger);

    if (TEfficiency::CheckConsistency(h_track_selection_cut, h_trigger))
        tr_eff_track_selection_cut = std::make_shared<TEfficiency>(h_track_selection_cut, h_trigger);

    if (TEfficiency::CheckConsistency(h_xtrl_cut, h_trigger))
        tr_eff_xtrl_cut = std::make_shared<TEfficiency>(h_xtrl_cut, h_trigger);

    if (TEfficiency::CheckConsistency(h_psd_charge_cut, h_trigger))
        tr_eff_psd_charge_cut = std::make_shared<TEfficiency>(h_psd_charge_cut, h_trigger);

    if (TEfficiency::CheckConsistency(h_stk_charge_cut, h_trigger))
        tr_eff_stk_charge_cut = std::make_shared<TEfficiency>(h_stk_charge_cut, h_trigger);

    if (TEfficiency::CheckConsistency(h_all_cut, h_trigger))
        tr_eff_all_cut = std::make_shared<TEfficiency>(h_all_cut, h_trigger);

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
    tr_eff_stk_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
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
    tr_eff_stk_charge_cut->SetName("tr_eff_stk_charge_cut");
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
    tr_eff_stk_charge_cut->SetTitle("stk charge cut efficiency");
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
    tr_eff_stk_charge_cut->Write();
    tr_eff_all_cut->Write();

    // Create geometric folder
    auto geometric_dir = ratioDir->mkdir("Geometric");
    geometric_dir->cd();

    // Define TEfficiency pointers
    std::shared_ptr<TEfficiency> geo_eff_maxElayer_cut;
    std::shared_ptr<TEfficiency> geo_eff_maxBarLayer_cut;
    std::shared_ptr<TEfficiency> geo_eff_BGOTrackContainment_cut;
    std::shared_ptr<TEfficiency> geo_eff_BGO_fiducial;
    std::shared_ptr<TEfficiency> geo_eff_nBarLayer13_cut;
    std::shared_ptr<TEfficiency> geo_eff_maxRms_cut;
    std::shared_ptr<TEfficiency> geo_eff_track_selection_cut;
    std::shared_ptr<TEfficiency> geo_eff_xtrl_cut;
    std::shared_ptr<TEfficiency> geo_eff_psd_charge_cut;
    std::shared_ptr<TEfficiency> geo_eff_stk_charge_cut;
    std::shared_ptr<TEfficiency> geo_eff_all_cut;

    if (TEfficiency::CheckConsistency(h_geometric_maxElayer_cut, h_gometric_cut))
        geo_eff_maxElayer_cut = std::make_shared<TEfficiency>(h_geometric_maxElayer_cut, h_gometric_cut);
    
    if (TEfficiency::CheckConsistency(h_geometric_maxBarLayer_cut, h_gometric_cut))
        geo_eff_maxBarLayer_cut = std::make_shared<TEfficiency>(h_geometric_maxBarLayer_cut, h_gometric_cut);

    if (TEfficiency::CheckConsistency(h_geometric_BGOTrackContainment_cut, h_gometric_cut))
        geo_eff_BGOTrackContainment_cut = std::make_shared<TEfficiency>(h_geometric_BGOTrackContainment_cut, h_gometric_cut);
    
    if (TEfficiency::CheckConsistency(h_geometric_BGO_fiducial_cut, h_gometric_cut))
        geo_eff_BGO_fiducial = std::make_shared<TEfficiency>(h_geometric_BGO_fiducial_cut, h_gometric_cut);

    if (TEfficiency::CheckConsistency(h_geometric_nBarLayer13_cut, h_gometric_cut))
        geo_eff_nBarLayer13_cut = std::make_shared<TEfficiency>(h_geometric_nBarLayer13_cut, h_gometric_cut);

    if (TEfficiency::CheckConsistency(h_geometric_maxRms_cut, h_gometric_cut))
        geo_eff_maxRms_cut = std::make_shared<TEfficiency>(h_geometric_maxRms_cut, h_gometric_cut);
    
    if (TEfficiency::CheckConsistency(h_geometric_track_selection_cut, h_gometric_cut))
        geo_eff_track_selection_cut = std::make_shared<TEfficiency>(h_geometric_track_selection_cut, h_gometric_cut);
    
    if (TEfficiency::CheckConsistency(h_geometric_xtrl_cut, h_gometric_cut))
        geo_eff_xtrl_cut = std::make_shared<TEfficiency>(h_geometric_xtrl_cut, h_gometric_cut);

    if (TEfficiency::CheckConsistency(h_geometric_psd_charge_cut, h_gometric_cut))
        geo_eff_psd_charge_cut = std::make_shared<TEfficiency>(h_geometric_psd_charge_cut, h_gometric_cut);

    if (TEfficiency::CheckConsistency(h_geometric_stk_charge_cut, h_gometric_cut))
        geo_eff_stk_charge_cut = std::make_shared<TEfficiency>(h_geometric_stk_charge_cut, h_gometric_cut);
    
    if (TEfficiency::CheckConsistency(h_geometric_all_cut, h_gometric_cut))
        geo_eff_all_cut = std::make_shared<TEfficiency>(h_geometric_all_cut, h_gometric_cut);
    
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
    geo_eff_stk_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
    geo_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);

    geo_eff_maxElayer_cut->SetName("geo_eff_maxElayer_cut");
    geo_eff_maxBarLayer_cut->SetName("geo_eff_maxBarLayer_cut");
    geo_eff_BGOTrackContainment_cut->SetName("geo_eff_BGOTrackContainment_cut");
    geo_eff_BGO_fiducial->SetName("geo_eff_BGO_fiducial");
    geo_eff_nBarLayer13_cut->SetName("geo_eff_nBarLayer13_cut");
    geo_eff_maxRms_cut->SetName("geo_eff_maxRms_cut");
    geo_eff_track_selection_cut->SetName("geo_eff_track_selection_cut");
    geo_eff_xtrl_cut->SetName("geo_eff_xtrl_cut");
    geo_eff_psd_charge_cut->SetName("geo_eff_psd_charge_cut");
    geo_eff_stk_charge_cut->SetName("geo_eff_stk_charge_cut");
    geo_eff_all_cut->SetName("geo_eff_all_cut");

    geo_eff_maxElayer_cut->SetTitle("geometic maxElayer cut efficiency");
    geo_eff_maxBarLayer_cut->SetTitle("geometric maxBarLayer cut efficiency");
    geo_eff_BGOTrackContainment_cut->SetTitle("geometric BGOTrackContainment cut efficiency");
    geo_eff_BGO_fiducial->SetTitle("geometric BGO fiducial cut efficiency");
    geo_eff_nBarLayer13_cut->SetTitle("geometric nBarLayer13 cut efficiency");
    geo_eff_maxRms_cut->SetTitle("geometric maxRms cut efficiency");
    geo_eff_track_selection_cut->SetTitle("geometric track selection cut efficiency");
    geo_eff_xtrl_cut->SetTitle("geometric xtrl cut efficiency");
    geo_eff_psd_charge_cut->SetTitle("geometric psd charge cut efficiency");
    geo_eff_stk_charge_cut->SetTitle("geometric stk charge cut efficiency");
    geo_eff_all_cut->SetTitle("geometric all cut efficiency");

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
    geo_eff_stk_charge_cut->Write();
    geo_eff_all_cut->Write();

    // Create BGO_fiducial_volume folder
    auto BGOfiducial_dir = ratioDir->mkdir("BGO_fiducial_volume");
    BGOfiducial_dir->cd();

    // Define TEfficiency pointers
    std::shared_ptr<TEfficiency> BGOfiducial_eff_nBarLayer13_cut;
    std::shared_ptr<TEfficiency> BGOfiducial_eff_maxRms_cut;
    std::shared_ptr<TEfficiency> BGOfiducial_eff_track_selection_cut;
    std::shared_ptr<TEfficiency> BGOfiducial_eff_xtrl_cut;
    std::shared_ptr<TEfficiency> BGOfiducial_eff_psd_charge_cut;
    std::shared_ptr<TEfficiency> BGOfiducial_eff_stk_charge_cut;
    std::shared_ptr<TEfficiency> BGOfiducial_eff_all_cut;

    if (TEfficiency::CheckConsistency(h_BGOfiducial_nBarLayer13_cut, h_BGO_fiducial_cut))
        BGOfiducial_eff_nBarLayer13_cut = std::make_shared<TEfficiency>(h_BGOfiducial_nBarLayer13_cut, h_BGO_fiducial_cut);

    if (TEfficiency::CheckConsistency(h_BGOfiducial_maxRms_cut, h_BGO_fiducial_cut))
        BGOfiducial_eff_maxRms_cut = std::make_shared<TEfficiency>(h_BGOfiducial_maxRms_cut, h_BGO_fiducial_cut);

    if (TEfficiency::CheckConsistency(h_BGOfiducial_track_selection_cut, h_BGO_fiducial_cut))
        BGOfiducial_eff_track_selection_cut = std::make_shared<TEfficiency>(h_BGOfiducial_track_selection_cut, h_BGO_fiducial_cut);

    if (TEfficiency::CheckConsistency(h_BGOfiducial_xtrl_cut, h_BGO_fiducial_cut))
        BGOfiducial_eff_xtrl_cut = std::make_shared<TEfficiency>(h_BGOfiducial_xtrl_cut, h_BGO_fiducial_cut);

    if (TEfficiency::CheckConsistency(h_BGOfiducial_psd_charge_cut, h_BGO_fiducial_cut))
        BGOfiducial_eff_psd_charge_cut = std::make_shared<TEfficiency>(h_BGOfiducial_psd_charge_cut, h_BGO_fiducial_cut);

    if (TEfficiency::CheckConsistency(h_BGOfiducial_stk_charge_cut, h_BGO_fiducial_cut))
        BGOfiducial_eff_stk_charge_cut = std::make_shared<TEfficiency>(h_BGOfiducial_stk_charge_cut, h_BGO_fiducial_cut);

    if (TEfficiency::CheckConsistency(h_BGOfiducial_all_cut, h_BGO_fiducial_cut))
        BGOfiducial_eff_all_cut = std::make_shared<TEfficiency>(h_BGOfiducial_all_cut, h_BGO_fiducial_cut);
    
    // Set uniform statistic option
    BGOfiducial_eff_nBarLayer13_cut->SetStatisticOption(TEfficiency::kBUniform);
    BGOfiducial_eff_maxRms_cut->SetStatisticOption(TEfficiency::kBUniform);
    BGOfiducial_eff_track_selection_cut->SetStatisticOption(TEfficiency::kBUniform);
    BGOfiducial_eff_xtrl_cut->SetStatisticOption(TEfficiency::kBUniform);
    BGOfiducial_eff_psd_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
    BGOfiducial_eff_stk_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
    BGOfiducial_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);

    BGOfiducial_eff_nBarLayer13_cut->SetName("BGOfiducial_eff_nBarLayer13_cut");
    BGOfiducial_eff_maxRms_cut->SetName("BGOfiducial_eff_maxRms_cut");
    BGOfiducial_eff_track_selection_cut->SetName("BGOfiducial_eff_track_selection_cut");
    BGOfiducial_eff_xtrl_cut->SetName("BGOfiducial_eff_xtrl_cut");
    BGOfiducial_eff_psd_charge_cut->SetName("BGOfiducial_eff_psd_charge_cut");
    BGOfiducial_eff_stk_charge_cut->SetName("BGOfiducial_eff_stk_charge_cut");
    BGOfiducial_eff_all_cut->SetName("BGOfiducial_eff_all_cut");

    BGOfiducial_eff_nBarLayer13_cut->SetName("BGOfiducial nBarLayer13 cut efficiency");
    BGOfiducial_eff_maxRms_cut->SetName("BGOfiducial maxRms cut efficiency");
    BGOfiducial_eff_track_selection_cut->SetName("BGOfiducial track selection cut efficiency");
    BGOfiducial_eff_xtrl_cut->SetName("BGOfiducial xtrl cut efficiency");
    BGOfiducial_eff_psd_charge_cut->SetName("BGOfiducial psd charge cut efficiency");
    BGOfiducial_eff_stk_charge_cut->SetName("BGOfiducial stk charge cut efficiency");
    BGOfiducial_eff_all_cut->SetName("BGOfiducial all cut efficiency");

    // Write histos to disk
    BGOfiducial_eff_nBarLayer13_cut->Write();
    BGOfiducial_eff_maxRms_cut->Write();
    BGOfiducial_eff_track_selection_cut->Write();
    BGOfiducial_eff_xtrl_cut->Write();
    BGOfiducial_eff_psd_charge_cut->Write();
    BGOfiducial_eff_stk_charge_cut->Write();
    BGOfiducial_eff_all_cut->Write();

    // Create output analysis dir in the output TFile
    auto preGeo_analysisDir = outFile.mkdir("Analysis_preGeoCut");
    preGeo_analysisDir->cd();

    h_preGeo_BGOrec_topX_vs_realX.Write();
    h_preGeo_BGOrec_topY_vs_realY.Write();
    h_preGeo_real_slopeX.Write();
    h_preGeo_real_slopeY.Write();
    h_preGeo_BGOrec_slopeX.Write();
    h_preGeo_BGOrec_slopeY.Write();
    h_preGeo_real_interceptX.Write();
    h_preGeo_real_interceptY.Write();
    h_preGeo_BGOrec_interceptX.Write();
    h_preGeo_BGOrec_interceptY.Write();
    h_preGeo_real_topMap.Write();
    h_preGeo_BGOreco_topMap.Write();
    h_preGeo_real_bottomMap.Write();
    h_preGeo_BGOreco_bottomMap.Write();

    h_noBGOenergy_real_topMap.Write();

    auto geo_analysisDir = outFile.mkdir("Analysis_GeoCut");
    geo_analysisDir->cd();

    h_geo_BGOrec_topX_vs_realX.Write();
    h_geo_BGOrec_topY_vs_realY.Write();
    h_geo_real_slopeX.Write();
    h_geo_real_slopeY.Write();
    h_geo_BGOrec_slopeX.Write();
    h_geo_BGOrec_slopeY.Write();
    h_geo_real_interceptX.Write();
    h_geo_real_interceptY.Write();
    h_geo_BGOrec_interceptX.Write();
    h_geo_BGOrec_interceptY.Write();
    h_geo_real_topMap.Write();
    h_geo_BGOreco_topMap.Write();
    h_geo_real_bottomMap.Write();
    h_geo_BGOreco_bottomMap.Write();

    auto BGOdir = outFile.mkdir("BGO_Energy");
    BGOdir->cd();

    h_BGOrec_E.Write();
    h_BGOrec_E_corr.Write();
    h_simu_energy.Write();
    h_energy_diff.Write();
    h_layer_max_energy_ratio.Write();

    h_triggered_BGOrec_E.Write();
    h_triggered_BGOrec_E_corr.Write();
    h_triggered_simu_energy.Write();
    h_triggered_energy_diff.Write();

    for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
        h_layer_energy_ratio[lIdx].Write();

    auto XTRLdir = outFile.mkdir("xtrl");
    XTRLdir->cd();

    h_xtrl_energy_int.Write();
    h_xtrl.Write();

    auto chargeDir = outFile.mkdir("STKcharge");
    chargeDir->cd();

    h_chargeX.Write();
    h_chargeY.Write();
    h_charge.Write();
    h_charge2D.Write();

    h_selected_chargeX.Write();
    h_selected_chargeY.Write();
    h_selected_charge.Write();
    h_selected_charge2D.Write();

    // Create ancillary output folder
    if (ancillary_cuts.compute_proton_background)
    {
        auto ancillaryDir = outFile.mkdir("mc_ancillary");
        ancillaryDir->cd();

        h_background_under_xtrl_cut.Write();
        h_background_over_xtrl_cut.Write();
        
        // Create proton background ratio
        auto proton_background_ratio = static_cast<TH1D*>(h_background_under_xtrl_cut.Clone("proton_background_ratio"));
        proton_background_ratio->SetTitle("Proton background ratio");
        proton_background_ratio->Divide(&h_background_over_xtrl_cut);

        proton_background_ratio->Write();
    }
    
    outFile.cd();
}