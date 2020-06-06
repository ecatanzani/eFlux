#include "acceptance.h"
#include "acceptance_cuts.h"
#include "data_cuts.h"
#include "energy_match.h"
#include "aggregate_events.h"
#include "wtsydp.h"
#include "working_dir.h"
#include "BGO_energy_cuts.h"
#include "DAMPE_geo_structure.h"
#include "DmpBgoContainer.h"
#include "read_sets_config_file.h"

#include "TGraphAsymmErrors.h"

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
    TH1D h_all_cut("h_all_cut", "Energy Distribution - All cut ", logEBins.size() - 1, &(logEBins[0]));

    // Acceptance - Cuts && Geometric Cut
    TH1D h_geometric_maxElayer_cut("h_geometric_maxElayer_cut", "Energy Distribution - maxElayer + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_maxBarLayer_cut("h_geometric_maxBarLayer_cut", "Energy Distribution - maxBarLayer + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_BGOTrackContainment_cut("h_geometric_BGOTrackContainment_cut", "Energy Distribution - BGOTrackContainment + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_BGO_fiducial("h_geometric_BGO_fiducial", "Energy Distibution - BGO fiducial + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_nBarLayer13_cut("h_geometric_nBarLayer13_cut", "Energy Distribution - nBarLayer13 + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_maxRms_cut("h_geometric_maxRms_cut", "Energy Distribution - maxRms + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_track_selection_cut("h_geometric_track_selection_cut", "Energy Distribution - track selection + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_xtrl_cut("h_geometric_xtrl_cut", "Energy Distribution - xtrl + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_psd_charge_cut("h_geometric_psd_charge_cut", "Energy Distribution - psd charge + geometric (trigger selection) cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_geometric_all_cut("h_geometric_all_cut", "Energy Distribution - All + geometric (trigger selection) cut ", logEBins.size() - 1, &(logEBins[0]));

    // Acceptance - Cuts && BGO fiducial volume cut
    TH1D h_BGOfiducial_nBarLayer13_cut("h_BGOfiducial_nBarLayer13_cut", "Energy Distribution - nBarLayer13 + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_maxRms_cut("h_BGOfiducial_maxRms_cut", "Energy Distribution - maxRms  + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_track_selection_cut("h_BGOfiducial_track_selection_cut", "Energy Distribution - track selection + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_xtrl_cut("h_BGOfiducial_xtrl_cut", "Energy Distribution - xtrl + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_psd_charge_cut("h_BGOfiducial_psd_charge_cut", "Energy Distribution - psd charge + BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOfiducial_all_cut("h_BGOfiducial_all_cut", "Energy Distribution - All + BGO fiducial cut ", logEBins.size() - 1, &(logEBins[0]));

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
    h_all_cut.Sumw2();

    // Sumw2 Acceptance - Cuts && Geometric Cut
    h_geometric_maxElayer_cut.Sumw2();
    h_geometric_maxBarLayer_cut.Sumw2();
    h_geometric_BGOTrackContainment_cut.Sumw2();
    h_geometric_BGO_fiducial.Sumw2();
    h_geometric_nBarLayer13_cut.Sumw2();
    h_geometric_maxRms_cut.Sumw2();
    h_geometric_track_selection_cut.Sumw2();
    h_geometric_xtrl_cut.Sumw2();
    h_geometric_psd_charge_cut.Sumw2();
    h_geometric_all_cut.Sumw2();

    // Sumw2 Acceptance - Cuts && BGO fiducial volume cut
    h_BGOfiducial_nBarLayer13_cut.Sumw2();
    h_BGOfiducial_maxRms_cut.Sumw2();
    h_BGOfiducial_track_selection_cut.Sumw2();
    h_BGOfiducial_xtrl_cut.Sumw2();
    h_BGOfiducial_psd_charge_cut.Sumw2();
    h_BGOfiducial_all_cut.Sumw2();

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

    // Create and load acceptance events cuts from config file
    cuts_conf acceptance_cuts;
    data_active_cuts active_cuts;
    load_acceptance_struct(acceptance_cuts, active_cuts, wd);

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

        // Fill the energy histos
        h_BGOrec_E.Fill(bgoTotalE_raw * _GeV);
        h_BGOrec_E_corr.Fill(bgoTotalE * _GeV);
        h_simu_energy.Fill(simuEnergy * _GeV);
        h_energy_diff.Fill((simuEnergy - bgoTotalE) * _GeV);

        // Don't accept events outside the selected energy window
        if (simuEnergy * _GeV < acceptance_cuts.min_event_energy || simuEnergy * _GeV > acceptance_cuts.max_event_energy)
            continue;

        if(geometric_cut(simu_primaries))
            h_geo_factor.Fill(simuEnergy * _GeV);

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
                h_incoming.Fill(simuEnergy * _GeV);
            }
            else
                continue;
        }
        else
        {
            // Increase the incoming events
            h_incoming.Fill(simuEnergy * _GeV);

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
                acceptance_cuts);
            filter_all_cut *= filter_xtrl_cut;
        }

        // **** psd_charge cut ****
        if (active_cuts.psd_charge)
        {
            filter_psd_charge_cut = psd_charge_cut(
                psdhits,
                bgorec,
                acceptance_cuts,
                event_best_track);
            filter_all_cut *= filter_psd_charge_cut;
        }

        // Fill cuts histos

        // Fill geometric cut histos
        if (filter_geometric_cut)
        {
            h_gometric_cut.Fill(simuEnergy * _GeV);

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
                h_geometric_maxElayer_cut.Fill(simuEnergy * _GeV);

            // Geometric cut && maxBarLayer cut
            if (filter_BGO_fiducial_maxBarLayer_cut)
                h_geometric_maxBarLayer_cut.Fill(simuEnergy * _GeV);
            
            // Geometric cut && BGOTrackContainment cut
            if (filter_BGO_fiducial_BGOTrackContainment_cut)
                h_geometric_BGOTrackContainment_cut.Fill(simuEnergy * _GeV);

            // Geometric cut && BGO fiducial cut
            if (filter_BGO_fiducial_cut)
                h_geometric_BGO_fiducial.Fill(simuEnergy * _GeV);

            // Geometric cut && nBarLayer13 cut
            if (filter_nBarLayer13_cut)
                h_geometric_nBarLayer13_cut.Fill(simuEnergy * _GeV);

            // Geometric cut && maxRms cut
            if (filter_maxRms_cut)
                h_geometric_maxRms_cut.Fill(simuEnergy * _GeV);

            // Geometric cut && track selection cut
            if (filter_track_selection_cut)
                h_geometric_track_selection_cut.Fill(simuEnergy * _GeV);

            // Geometric cut && XTRL cut
            if (filter_xtrl_cut)
                h_geometric_xtrl_cut.Fill(simuEnergy * _GeV);

            // Geometric cut && PSD charge cut
            if (filter_psd_charge_cut)
                h_geometric_psd_charge_cut.Fill(simuEnergy * _GeV);

            // Geometric cut and all cuts
            if (filter_all_cut)
                h_geometric_all_cut.Fill(simuEnergy * _GeV);
        }

        // Fill BGO_fiducial_maxElayer cut histo
        if (filter_BGO_fiducial_maxElayer_cut)
            h_maxElayer_cut.Fill(simuEnergy * _GeV);

        // Fill BGO_fiducial_maxBarLayer cut histo
        if (filter_BGO_fiducial_maxBarLayer_cut)
            h_maxBarLayer_cut.Fill(simuEnergy * _GeV);

        // Fill BGO_fiducial_BGOTrackContainment cut histo
        if (filter_BGO_fiducial_BGOTrackContainment_cut)
            h_BGOTrackContainment_cut.Fill(simuEnergy * _GeV);

        // Fill BGO fiducial volume cut
        if (filter_BGO_fiducial_cut)
        {
            h_BGO_fiducial_cut.Fill(simuEnergy * _GeV);
        
            // BGO fiducial cut && nBarLayer13 cut
            if (filter_nBarLayer13_cut)
                h_BGOfiducial_nBarLayer13_cut.Fill(simuEnergy * _GeV);
            
            // BGO fiducial cut && maxRms cut
            if (filter_maxRms_cut)
                h_BGOfiducial_maxRms_cut.Fill(simuEnergy * _GeV);

            // BGO fiducial cut && track selection cut
            if (filter_track_selection_cut)
                h_BGOfiducial_track_selection_cut.Fill(simuEnergy * _GeV);

            // BGO fiducial cut && XTRL cut
            if (filter_xtrl_cut)
                h_BGOfiducial_xtrl_cut.Fill(simuEnergy * _GeV);

            // BGO fiducial cut && PSD charge cut
            if (filter_psd_charge_cut)
                h_BGOfiducial_psd_charge_cut.Fill(simuEnergy * _GeV);

            // BGO fiducial cut && all cut
            if (filter_all_cut)
                h_BGOfiducial_all_cut.Fill(simuEnergy * _GeV);
        }

        // Fill nBarLayer13 cut histo
        if (filter_nBarLayer13_cut)
            h_nBarLayer13_cut.Fill(simuEnergy * _GeV);

        // Fill maxRms cut histo
        if (filter_maxRms_cut)
            h_maxRms_cut.Fill(simuEnergy * _GeV);

        // Fill track selection cut histo
        if (filter_track_selection_cut)
            h_track_selection_cut.Fill(simuEnergy * _GeV);

        // Fill XTRL cut histo
        if (filter_xtrl_cut)
            h_xtrl_cut.Fill(simuEnergy * _GeV);

        // Fill PSD charge cut histo
        if (filter_psd_charge_cut)
            h_psd_charge_cut.Fill(simuEnergy * _GeV);

        // Fill all cut histo
        if (active_cuts.nActiveCuts)
            if (filter_all_cut)
                h_all_cut.Fill(simuEnergy * _GeV);
    }

    if (verbose)
    {
        std::cout << "\n\n ****** \n\n";
        std::cout << "generated events in good energy range: " << h_incoming.GetEntries() << std::endl;
        std::cout << "triggered events: " << h_trigger.GetEntries() << std::endl;

        auto refEntries = h_trigger.GetEntries();

        if(h_geo_factor.GetEntries())
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
    std::vector<double> acceptanceError_all_cut(acceptanceError_geometric_factor.size(), 0);
    
    std::vector<double> energy_LowError(energyValues.size(), 0);
    std::vector<double> energy_HighError(energyValues.size(), 0);

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

        acceptanceError_geometric_factor[index] = h_acceptance_geometric_factor->GetBinError(index+1)/2.;
        acceptanceError_gometric_cut[index] = h_acceptance_gometric_cut->GetBinError(index+1)/2.;
        acceptanceError_maxElayer_cut[index] = h_acceptance_maxElayer_cut->GetBinError(index+1)/2.;
        acceptanceError_maxBarLayer_cut[index] = h_acceptance_maxBarLayer_cut->GetBinError(index+1)/2.;
        acceptanceError_BGOTrackContainment_cut[index] = h_acceptance_BGOTrackContainment_cut->GetBinError(index+1)/2.;
        acceptanceError_BGO_fiducial_cut[index] = h_acceptance_BGO_fiducial_cut->GetBinError(index+1)/2.;
        acceptanceError_BGO_fiducial_nBarLayer13_cut[index] = h_acceptance_BGO_fiducial_nBarLayer13_cut->GetBinError(index+1)/2.;
        acceptanceError_BGO_fiducial_maxRms_cut[index] = h_acceptance_BGO_fiducial_maxRms_cut->GetBinError(index+1)/2.;
        acceptanceError_BGO_fiducial_track_selection_cut[index] = h_acceptance_BGO_fiducial_track_selection_cut->GetBinError(index+1)/2.;
        acceptanceError_BGO_fiducial_xtrl_cut[index] = h_acceptance_BGO_fiducial_xtrl_cut->GetBinError(index+1)/2.;
        acceptanceError_BGO_fiducial_psd_charge_cut[index] = h_acceptance_BGO_fiducial_psd_charge_cut->GetBinError(index+1)/2.;
        acceptanceError_nBarLayer13_cut[index] = h_acceptance_nBarLayer13_cut->GetBinError(index+1)/2.;
        acceptanceError_maxRms_cut[index] = h_acceptance_maxRms_cut->GetBinError(index+1)/2.;
        acceptanceError_track_selection_cut[index] = h_acceptance_track_selection_cut->GetBinError(index+1)/2.;
        acceptanceError_xtrl_cut[index] = h_acceptance_xtrl_cut->GetBinError(index+1)/2.;
        acceptanceError_psd_charge_cut[index] = h_acceptance_psd_charge_cut->GetBinError(index+1)/2.;
        acceptanceError_all_cut[index] = h_acceptance_all_cut->GetBinError(index+1)/2.; 

        energy_LowError[index] = energyValues[index] - *it;
        energy_HighError[index] = *(it + 1) - energyValues[index];
    }

    // Building graphs
     TGraphAsymmErrors gr_acceptance_geometric_factor(energyValues.size(), &energyValues[0], &acceptanceValues_geometric_factor[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_geometric_factor[0], &acceptanceError_geometric_factor[0]);
    TGraphAsymmErrors gr_acceptance_gometric_cut(energyValues.size(), &energyValues[0], &acceptanceValues_gometric_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_gometric_cut[0], &acceptanceError_gometric_cut[0]);
    TGraphAsymmErrors gr_acceptance_maxElayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxElayer_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_maxElayer_cut[0], &acceptanceError_maxElayer_cut[0]);
    TGraphAsymmErrors gr_acceptance_maxBarLayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxBarLayer_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_maxBarLayer_cut[0], &acceptanceError_maxBarLayer_cut[0]);
    TGraphAsymmErrors gr_acceptance_BGOTrackContainment_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGOTrackContainment_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_BGOTrackContainment_cut[0], &acceptanceError_BGOTrackContainment_cut[0]);
    TGraphAsymmErrors gr_acceptance_BGO_fiducial_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_BGO_fiducial_cut[0], &acceptanceError_BGO_fiducial_cut[0]);
    TGraphAsymmErrors gr_acceptance_BGO_fiducial_nBarLayer13_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_nBarLayer13_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_BGO_fiducial_nBarLayer13_cut[0], &acceptanceError_BGO_fiducial_nBarLayer13_cut[0]);
    TGraphAsymmErrors gr_acceptance_BGO_fiducial_maxRms_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_maxRms_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_BGO_fiducial_maxRms_cut[0], &acceptanceError_BGO_fiducial_maxRms_cut[0]);
    TGraphAsymmErrors gr_acceptance_BGO_fiducial_track_selection_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_track_selection_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_BGO_fiducial_track_selection_cut[0], &acceptanceError_BGO_fiducial_track_selection_cut[0]);
    TGraphAsymmErrors gr_acceptance_BGO_fiducial_xtrl_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_xtrl_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_BGO_fiducial_xtrl_cut[0], &acceptanceError_BGO_fiducial_xtrl_cut[0]);
    TGraphAsymmErrors gr_acceptance_BGO_fiducial_psd_charge_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_psd_charge_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_BGO_fiducial_psd_charge_cut[0], &acceptanceError_BGO_fiducial_psd_charge_cut[0]);
    TGraphAsymmErrors gr_acceptance_nBarLayer13_cut(energyValues.size(), &energyValues[0], &acceptanceValues_nBarLayer13_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_nBarLayer13_cut[0], &acceptanceError_nBarLayer13_cut[0]);
    TGraphAsymmErrors gr_acceptance_maxRms_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxRms_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_maxRms_cut[0], &acceptanceError_maxRms_cut[0]);
    TGraphAsymmErrors gr_acceptance_track_selection_cut(energyValues.size(), &energyValues[0], &acceptanceValues_track_selection_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_track_selection_cut[0], &acceptanceError_track_selection_cut[0]);
    TGraphAsymmErrors gr_acceptance_xtrl_cut(energyValues.size(), &energyValues[0], &acceptanceValues_xtrl_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_xtrl_cut[0], &acceptanceError_xtrl_cut[0]);
    TGraphAsymmErrors gr_acceptance_psd_charge_cut(energyValues.size(), &energyValues[0], &acceptanceValues_psd_charge_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_psd_charge_cut[0], &acceptanceError_psd_charge_cut[0]);
    TGraphAsymmErrors gr_acceptance_all_cut(energyValues.size(), &energyValues[0], &acceptanceValues_all_cut[0], &energy_LowError[0], &energy_HighError[0], &acceptanceError_all_cut[0], &acceptanceError_all_cut[0]);

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
    h_all_cut.Write();

    // Acceptance - Cuts && Geometric Cut
    h_geometric_maxElayer_cut.Write();
    h_geometric_maxBarLayer_cut.Write();
    h_geometric_BGOTrackContainment_cut.Write();
    h_geometric_BGO_fiducial.Write();
    h_geometric_nBarLayer13_cut.Write();
    h_geometric_maxRms_cut.Write();
    h_geometric_track_selection_cut.Write();
    h_geometric_xtrl_cut.Write();
    h_geometric_psd_charge_cut.Write();
    h_geometric_all_cut.Write();
    // Acceptance - Cuts && BGO fiducial volume cut
    h_BGOfiducial_nBarLayer13_cut.Write();
    h_BGOfiducial_maxRms_cut.Write();
    h_BGOfiducial_track_selection_cut.Write();
    h_BGOfiducial_xtrl_cut.Write();
    h_BGOfiducial_psd_charge_cut.Write();
    h_BGOfiducial_all_cut.Write();

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

    // Building ratio histos
    auto h_trigger_efficiency = static_cast<TH1D *>(h_gometric_cut.Clone("h_trigger_efficiency"));
    auto h_ratio_tr_gometric_cut = static_cast<TH1D *>(h_gometric_cut.Clone("h_ratio_tr_gometric_cut"));
    auto h_ratio_tr_maxElayer_cut = static_cast<TH1D *>(h_maxElayer_cut.Clone("h_ratio_tr_maxElayer_cut"));
    auto h_ratio_tr_maxBarLayer_cut = static_cast<TH1D *>(h_maxBarLayer_cut.Clone("h_ratio_tr_maxBarLayer_cut"));
    auto h_ratio_tr_BGOTrackContainment_cut = static_cast<TH1D *>(h_BGOTrackContainment_cut.Clone("h_ratio_tr_BGOTrackContainment_cut"));
    auto h_ratio_tr_BGO_fiducial = static_cast<TH1D *>(h_BGO_fiducial_cut.Clone("h_ratio_tr_BGO_fiducial"));
    auto h_ratio_tr_nBarLayer13_cut = static_cast<TH1D *>(h_nBarLayer13_cut.Clone("h_ratio_tr_nBarLayer13_cut"));
    auto h_ratio_tr_maxRms_cut = static_cast<TH1D *>(h_maxRms_cut.Clone("h_ratio_tr_maxRms_cut"));
    auto h_ratio_tr_track_selection_cut = static_cast<TH1D *>(h_track_selection_cut.Clone("h_ratio_tr_track_selection_cut"));
    auto h_ratio_tr_xtrl_cut = static_cast<TH1D *>(h_xtrl_cut.Clone("h_ratio_tr_xtrl_cut"));
    auto h_ratio_tr_psd_charge_cut = static_cast<TH1D *>(h_psd_charge_cut.Clone("h_ratio_tr_psd_charge_cut"));
    auto h_ratio_tr_all_cut = static_cast<TH1D *>(h_all_cut.Clone("h_ratio_tr_all_cut"));

    // Scale histos respect to the trigger cut events
    h_trigger_efficiency->Divide(&h_geo_factor);
    h_ratio_tr_gometric_cut->Divide(&h_trigger);
    h_ratio_tr_maxElayer_cut->Divide(&h_trigger);
    h_ratio_tr_maxBarLayer_cut->Divide(&h_trigger);
    h_ratio_tr_BGOTrackContainment_cut->Divide(&h_trigger);
    h_ratio_tr_BGO_fiducial->Divide(&h_trigger);
    h_ratio_tr_nBarLayer13_cut->Divide(&h_trigger);
    h_ratio_tr_maxRms_cut->Divide(&h_trigger);
    h_ratio_tr_track_selection_cut->Divide(&h_trigger);
    h_ratio_tr_xtrl_cut->Divide(&h_trigger);
    h_ratio_tr_psd_charge_cut->Divide(&h_trigger);
    h_ratio_tr_all_cut->Divide(&h_trigger);

    //Write histos to disk
    h_trigger_efficiency->Write();
    h_ratio_tr_gometric_cut->Write();
    h_ratio_tr_maxElayer_cut->Write();
    h_ratio_tr_maxBarLayer_cut->Write();
    h_ratio_tr_BGOTrackContainment_cut->Write();
    h_ratio_tr_BGO_fiducial->Write();
    h_ratio_tr_nBarLayer13_cut->Write();
    h_ratio_tr_maxRms_cut->Write();
    h_ratio_tr_track_selection_cut->Write();
    h_ratio_tr_xtrl_cut->Write();
    h_ratio_tr_psd_charge_cut->Write();
    h_ratio_tr_all_cut->Write();

    // Return to main ratio dir
    ratioDir->cd();

    // Create geometric folder
    auto geometric_dir = ratioDir->mkdir("Geometric");
    geometric_dir->cd();

    // Building ratio histos
    auto h_ratio_geo_maxElayer_cut = static_cast<TH1D *>(h_geometric_maxElayer_cut.Clone("h_ratio_geo_maxElayer_cut"));
    auto h_ratio_geo_maxBarLayer_cut = static_cast<TH1D *>(h_geometric_maxBarLayer_cut.Clone("h_ratio_geo_maxBarLayer_cut"));
    auto h_ratio_geo_BGOTrackContainment_cut = static_cast<TH1D *>(h_geometric_BGOTrackContainment_cut.Clone("h_ratio_geo_BGOTrackContainment_cut"));
    auto h_ratio_geo_BGO_fiducial = static_cast<TH1D *>(h_geometric_BGO_fiducial.Clone("h_ratio_geo_BGO_fiducial"));
    auto h_ratio_geo_nBarLayer13_cut = static_cast<TH1D *>(h_geometric_nBarLayer13_cut.Clone("h_ratio_geo_nBarLayer13_cut"));
    auto h_ratio_geo_maxRms_cut = static_cast<TH1D *>(h_geometric_maxRms_cut.Clone("h_ratio_geo_maxRms_cut"));
    auto h_ratio_geo_track_selection_cut = static_cast<TH1D *>(h_geometric_track_selection_cut.Clone("h_ratio_geo_track_selection_cut"));
    auto h_ratio_geo_xtrl_cut = static_cast<TH1D *>(h_geometric_xtrl_cut.Clone("h_ratio_geo_xtrl_cut"));
    auto h_ratio_geo_psd_charge_cut = static_cast<TH1D *>(h_geometric_psd_charge_cut.Clone("h_ratio_geo_psd_charge_cut"));
    auto h_ratio_geo_all_cut = static_cast<TH1D *>(h_geometric_all_cut.Clone("h_ratio_geo_all_cut"));

    // Scale histos respect to the geometric cut events
    h_ratio_geo_maxElayer_cut->Divide(&h_gometric_cut);
    h_ratio_geo_maxBarLayer_cut->Divide(&h_gometric_cut);
    h_ratio_geo_BGOTrackContainment_cut->Divide(&h_gometric_cut);
    h_ratio_geo_BGO_fiducial->Divide(&h_gometric_cut);
    h_ratio_geo_nBarLayer13_cut->Divide(&h_gometric_cut);
    h_ratio_geo_maxRms_cut->Divide(&h_gometric_cut);
    h_ratio_geo_track_selection_cut->Divide(&h_gometric_cut);
    h_ratio_geo_xtrl_cut->Divide(&h_gometric_cut);
    h_ratio_geo_psd_charge_cut->Divide(&h_gometric_cut);
    h_ratio_geo_all_cut->Divide(&h_gometric_cut);

    //Write histos to disk
    h_ratio_geo_maxElayer_cut->Write();
    h_ratio_geo_maxBarLayer_cut->Write();
    h_ratio_geo_BGOTrackContainment_cut->Write();
    h_ratio_geo_BGO_fiducial->Write();
    h_ratio_geo_nBarLayer13_cut->Write();
    h_ratio_geo_maxRms_cut->Write();
    h_ratio_geo_track_selection_cut->Write();
    h_ratio_geo_xtrl_cut->Write();
    h_ratio_geo_psd_charge_cut->Write();
    h_ratio_geo_all_cut->Write();

    auto BGOfiducial_dir = ratioDir->mkdir("BGO_fiducial_volume");
    BGOfiducial_dir->cd();

    // Building ratio histos
    auto h_ratio_BGOfiducial_nBarLayer13_cut = static_cast<TH1D *>(h_BGOfiducial_nBarLayer13_cut.Clone("h_ratio_BGOfiducial_nBarLayer13_cut"));
    auto h_ratio_BGOfiducial_maxRms_cut = static_cast<TH1D *>(h_BGOfiducial_maxRms_cut.Clone("h_ratio_BGOfiducial_maxRms_cut"));
    auto h_ratio_BGOfiducial_track_selection_cut = static_cast<TH1D *>(h_BGOfiducial_track_selection_cut.Clone("h_ratio_BGOfiducial_track_selection_cut"));
    auto h_ratio_BGOfiducial_xtrl_cut = static_cast<TH1D *>(h_BGOfiducial_xtrl_cut.Clone("h_ratio_BGOfiducial_xtrl_cut"));
    auto h_ratio_BGOfiducial_psd_charge_cut = static_cast<TH1D *>(h_BGOfiducial_psd_charge_cut.Clone("h_ratio_BGOfiducial_psd_charge_cut"));
    auto h_ratio_BGOfiducial_all_cut = static_cast<TH1D *>(h_BGOfiducial_all_cut.Clone("h_ratio_BGOfiducial_all_cut"));

    // Scale histos respect to the BGO fiducial cut events
    h_ratio_BGOfiducial_nBarLayer13_cut->Divide(&h_BGO_fiducial_cut);
    h_ratio_BGOfiducial_maxRms_cut->Divide(&h_BGO_fiducial_cut);
    h_ratio_BGOfiducial_track_selection_cut->Divide(&h_BGO_fiducial_cut);
    h_ratio_BGOfiducial_xtrl_cut->Divide(&h_BGO_fiducial_cut);
    h_ratio_BGOfiducial_psd_charge_cut->Divide(&h_BGO_fiducial_cut);
    h_ratio_BGOfiducial_all_cut->Divide(&h_BGO_fiducial_cut);

    //Write histos to disk
    h_ratio_BGOfiducial_nBarLayer13_cut->Write();
    h_ratio_BGOfiducial_maxRms_cut->Write();
    h_ratio_BGOfiducial_track_selection_cut->Write();
    h_ratio_BGOfiducial_xtrl_cut->Write();
    h_ratio_BGOfiducial_psd_charge_cut->Write();
    h_ratio_BGOfiducial_all_cut->Write();

    outFile.cd();

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

    outFile.cd();

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

    outFile.cd();

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

    outFile.cd();
}

#if 0
void buildAcceptance_vector(
    const std::string accInputPath,
    const bool verbose,
    const std::vector<float> &logEBins,
    TFile &outFile,
    const std::string wd)
{

    //auto dmpch = aggregateEventsDmpChain(accInputPath,verbose);
    auto dmpch = aggregateEventsTChain(accInputPath, verbose);

    // SimuPrimaries container
    std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries = std::make_shared<DmpEvtSimuPrimaries>();
    dmpch->SetBranchAddress("DmpEvtSimuPrimaries", &simu_primaries);

    // Register BGO constainer
    std::shared_ptr<DmpEvtBgoHits> bgohits = std::make_shared<DmpEvtBgoHits>();
    dmpch->SetBranchAddress("DmpEvtBgoHits", &bgohits);

    // Register BGO REC constainer
    std::shared_ptr<DmpEvtBgoRec> bgorec = std::make_shared<DmpEvtBgoRec>();
    dmpch->SetBranchAddress("DmpEvtBgoRec", &bgorec);

    // Event loop
    auto nevents = dmpch->GetEntries();
    if (verbose)
        std::cout << "\n\nTotal number of events: " << nevents << "\n\n";

    // Initialize the particle counter for each energy bin
    std::vector<unsigned int> gFactorCounts(logEBins.size() - 1, 0);
    std::vector<unsigned int> gen_gFactorCounts(logEBins.size() - 1, 0);

    // Create and load acceptance events cuts from config file
    acceptance_conf acceptance_cuts;
    load_acceptance_struct(acceptance_cuts, wd);

    double _GeV = 0.001;

    for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx)
    {
        // Get chain event
        dmpch->GetEvent(evIdx);

        // Event printout
        if (verbose)
            if (((evIdx + 1) % _kStep) == 0)
                std::cout << "\nProcessed " << evIdx + 1 << " events / " << nevents;

        // GOOD EVENT variable
        bool passEvent = true;

        // Get event total energy
        //double bgoTotalE = bgorec->GetTotalEnergy(); // Energy in MeV - not corrected
        double bgoTotalE = bgorec->GetElectronEcor();    // Returns corrected energy assuming this was an electron (MeV)
        double simuEnergy = simu_primaries->pvpart_ekin; //Energy of simu primaries particle in MeV

        // Don't accept events outside the selected energy window
        if (simuEnergy * _GeV < acceptance_cuts.min_event_energy || simuEnergy * _GeV > acceptance_cuts.max_event_energy)
            continue;

        allocateParticleEnergy(gen_gFactorCounts, logEBins, simuEnergy);

        // Don't process events that didn't hit the detector - for this I need to use the reco energy
        if (bgoTotalE == 0)
            continue;

        std::vector<std::vector<short>> layerBarIndex(DAMPE_bgo_nLayers, std::vector<short>());  // arrange BGO hits by layer
        std::vector<std::vector<short>> layerBarNumber(DAMPE_bgo_nLayers, std::vector<short>()); // arrange BGO bars by layer

        // Get the number of BGO hits
        int nBgoHits = bgohits->GetHittedBarNumber();

        // Scan BGO hits
        for (int ihit = 0; ihit < nBgoHits; ++ihit)
        {
            // Get layer ID
            auto layerID = bgohits->GetLayerID(ihit);
            // Get bar global ID
            auto iBar = ((bgohits->fGlobalBarID)[ihit] >> 6) & 0x1f;
            layerBarIndex[layerID].push_back(ihit);
            layerBarNumber[layerID].push_back(iBar);
        }

        std::vector<double> rmsLayer(DAMPE_bgo_nLayers, 0);
        std::vector<double> fracLayer(DAMPE_bgo_nLayers, 0);
        std::vector<double> eLayer(DAMPE_bgo_nLayers, 0);
        std::vector<double> eCoreLayer(DAMPE_bgo_nLayers, 0);
        std::vector<double> eCoreCoord(DAMPE_bgo_nLayers, 0);
        double sumRms = 0;

        for (int lay = 0; lay < DAMPE_bgo_nLayers; ++lay)
        {
            // Setting default value for maximum bar index and energy for each layer
            int imax = -1;
            if (layerBarIndex[lay].size())
                imax = layerBarIndex[lay].at(0);
            double maxE = (bgohits->fEnergy)[0];

            // Find the maximum of the nergy release in a certain layer, together with the bar ID
            for (auto it = layerBarNumber[lay].begin(); it != layerBarNumber[lay].end(); ++it)
            {
                int ihit = *it;
                double hitE = (bgohits->fEnergy)[ihit];
                if (hitE > maxE)
                {
                    maxE = hitE;
                    imax = ihit;
                }
            }

            rmsLayer[lay] = 0;
            fracLayer[lay] = 0;
            eLayer[lay] = 0;

            if (maxE)
            {
                // Find the bar index regarding the maximum energy release in a certain layer
                auto iBarMax = ((bgohits->fGlobalBarID)[imax] >> 6) & 0x1f;
                // Register the maximum energy release of a layer
                eCoreLayer[lay] = maxE;
                // Find the coordinate (weighted by the nergy release) of the bar with the biggest energy release in a certain layer
                auto coordMax = lay % 2 ? bgohits->GetHitX(imax) : bgohits->GetHitY(imax);
                eCoreCoord[lay] = maxE * coordMax;
                // Consider the nearest bar respect to the max one in order to better interpolate the position
                if (iBarMax > 0 && iBarMax < 21)
                {
                    for (auto it = layerBarNumber[lay].begin(); it != layerBarNumber[lay].end(); ++it)
                    {
                        auto ihit = *it;
                        auto iBar = ((bgohits->fGlobalBarID)[ihit] >> 6) & 0x1f;
                        if (iBar - iBarMax == 1 || iBar - iBarMax == -1)
                        {
                            double hitE = (bgohits->fEnergy)[ihit];
                            double thisCoord = lay % 2 ? bgohits->GetHitX(ihit) : bgohits->GetHitY(ihit);
                            eCoreLayer[lay] += hitE;
                            eCoreCoord[lay] += hitE * thisCoord;
                        }
                    }
                }
                // Get the CoG coordinate of the max energy bar cluster
                eCoreCoord[lay] /= eCoreLayer[lay];
                // Get the energy RMS of a layer
                for (auto it = layerBarNumber[lay].begin(); it != layerBarNumber[lay].end(); ++it)
                {
                    auto ihit = *it;
                    auto hitE = (bgohits->fEnergy)[ihit];
                    auto thisCoord = lay % 2 ? bgohits->GetHitX(ihit) : bgohits->GetHitY(ihit);
                    eLayer[lay] += hitE;
                    rmsLayer[lay] += hitE * (thisCoord - eCoreCoord[lay]) * (thisCoord - eCoreCoord[lay]);
                }
                rmsLayer[lay] = sqrt(rmsLayer[lay] / eLayer[lay]);
                fracLayer[lay] = eLayer[lay] / bgoTotalE;
                if (layerBarNumber[lay].size() <= 1)
                    rmsLayer[lay] = 0;
            }
            sumRms += rmsLayer[lay];
        }

        // Build XTR
        //auto Xtr = pow(sumRms, 4) * fracLayer[13] / 8000000.;

        /* ********************************* */

        if (maxElater_cut(bgorec, acceptance_cuts, bgoTotalE))
            if (maxBarLayer_cut(bgohits, nBgoHits))
                if (BGOTrackContainment_cut(bgorec, acceptance_cuts, passEvent))
                    // Use the simu energy as a reference to build the acceptance
                    allocateParticleEnergy(gFactorCounts, logEBins, simuEnergy);
        //std::cout << "\nSimuEnergy: " << simuEnergy << "\t RecoEnergy: " << bgoTotalE;
    }

    if (verbose)
    {
        unsigned int ev_counter = 0;
        std::cout << "\n\nFiltered events: " << std::accumulate(gFactorCounts.begin(), gFactorCounts.end(), ev_counter) << "/" << nevents << "\n\n";
    }
    for (auto it = gen_gFactorCounts.begin(); it != gen_gFactorCounts.end(); ++it)
        if (*it == 0)
        {
            auto index = std::distance(gen_gFactorCounts.begin(), it);
            gen_gFactorCounts[index] = 1;
        }
    double genSurface = 4 * TMath::Pi() * pow(acceptance_cuts.vertex_radius, 2) / 2;
    std::vector<double> d_gFactorCounts(gFactorCounts.begin(), gFactorCounts.end());
    std::vector<double> d_gen_gFactorCounts(gen_gFactorCounts.begin(), gen_gFactorCounts.end());
    std::transform(d_gFactorCounts.begin(), d_gFactorCounts.end(), d_gFactorCounts.begin(), [&genSurface](auto &elm) { return elm * genSurface; });
    std::vector<double> gFactor;
    const std::size_t vSize = std::min(d_gFactorCounts.size(), d_gen_gFactorCounts.size());
    std::transform(std::begin(d_gFactorCounts), std::begin(d_gFactorCounts) + vSize, std::begin(d_gen_gFactorCounts), std::back_inserter(gFactor), std::divides<double>{});
    std::vector<double> energyValues(gFactor.size(), 0);
    for (auto it = logEBins.begin(); it != (logEBins.end() - 1); ++it)
    {
        auto index = std::distance(logEBins.begin(), it);
        energyValues[index] = wtsydp(*it, *(it + 1), -1);
        //std::cout << std::endl << *it << "\t" << energyValues[index] << "\t" << *(it + 1) << std::endl;
    }
    if (!strcmp(_memType, "graph"))
    {
        // Create acceptance TGraph
        TGraph acceptanceGr(gFactor.size(), &energyValues[0], &gFactor[0]);
        acceptanceGr.SetName("acceptance");
        acceptanceGr.SetTitle("Acceptance");

        // Create output acceptance dir in the output TFile
        auto acceptanceDIr = outFile.mkdir("Acceptance");
        acceptanceDIr->cd();

        // Write final TGraph
        acceptanceGr.Write();
    }
    else if (!strcmp(_memType, "histo"))
    {
        auto acceptanceHisto = buildHistoFromVector(energyValues, gFactor);
        acceptanceHisto->SetName("acceptance");
        acceptanceHisto->SetTitle("acceptance");

        // Create output acceptance dir in the output TFile
        auto acceptanceDIr = outFile.mkdir("Acceptance");
        acceptanceDIr->cd();

        // Write final TGraph
        acceptanceHisto->Write();
    }
    else
    {
        std::cerr << "\n\nERROR: incorrect memtype: " << _memType;
        exit(123);
    }

    // Return to main TFile directory
    outFile.cd();
}
#endif