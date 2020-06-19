#include "data_loop.h"
#include "aggregate_events.h"
#include "read_sets_config_file.h"
#include "data_cuts.h"
#include "BGO_energy_cuts.h"
#include "flux.h"
#include "wtsydp.h"
#include "binning.h"

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

inline void updateProcessStatus(const int evIdx, int &kStep, const int nevents)
{
    auto percentage = ((evIdx + 1) / (double)nevents) * 100;
    if (floor(percentage) != 0 && ((int)floor(percentage) % kStep) == 0)
    {
        std::cout << "\n"
                  << floor(percentage) << " %\t | \tProcessed " << evIdx + 1 << " events / " << nevents;
        kStep += 10;
    }
}

TH1D evLoop(
    const std::vector<float> &logEBins,
    const std::string inputPath,
    TFile &outFile,
    const bool verbose,
    const std::string wd)
{
    //auto dmpch = aggregateEventsDmpChain(inputPath,verbose);
    auto dmpch = aggregateEventsTChain(inputPath, verbose);

    // Register Header container
    std::shared_ptr<DmpEvtHeader> evt_header = std::make_shared<DmpEvtHeader>();
    dmpch->SetBranchAddress("EventHeader", &evt_header);
    
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
    
    // Analysis histos - simu and reco energy of triggered events
    TH1D h_triggered_BGOrec_E("h_triggered_BGOrec_E", "Triggered BGO Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_triggered_BGOrec_E_corr("h_triggered_BGOrec_E_corr", "Triggered BGO Corrected Energy", logEBins.size() - 1, &(logEBins[0]));
    
    // After Geometric Cut
    // Slope X and Y
    TH1D h_geo_BGOrec_slopeX("h_geo_BGOrec_slopeX", "BGOrec Slope X", 1000, -90, 90);
    TH1D h_geo_BGOrec_slopeY("h_geo_BGOrec_slopeY", "BGOrec Slope Y", 1000, -90, 90);

    // Intercept X and Y
    TH1D h_geo_BGOrec_interceptX("h_geo_BGOrec_interceptX", "BGOrec Intercept X", 500, -500, 500);
    TH1D h_geo_BGOrec_interceptY("h_geo_BGOrec_interceptY", "BGOrec Intercept Y", 500, -500, 500);

    // Top Maps
    TH2D h_geo_BGOreco_topMap("h_geo_BGOreco_topMap", "BGOreco TOP Map", 500, -500, 500, 500, -500, 500);

    // Bottom Maps
    TH2D h_geo_BGOreco_bottomMap("h_geo_BGOreco_bottomMap", "BGOreco BOTTOM Map", 500, -500, 500, 500, -500, 500);

    // Ratio of layer energy respect to total BGO energy
    TH1D h_layer_max_energy_ratio("h_layer_max_energy_ratio", "Layer Energy Ratio", 100, 0, 1);
    std::vector<TH1D> h_layer_energy_ratio;
    init_BGO_histos(h_layer_energy_ratio);

    // XTRL histos
    auto xtrl_bins = LinearSpacedArray(0, 100, 1000);

    TH1D h_xtrl_energy_int("h_xtrl_energy_int","Energy integrated XTRL distribution", xtrl_bins.size() -1, &(xtrl_bins[0]));
    TH2D h_xtrl("h_xtrl", "XTRL energy Distribution", logEBins.size() - 1, &(logEBins[0]), xtrl_bins.size(), &(xtrl_bins[0]));

    // Sumw2 Acceptance - First-Cut histos
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
    
    // Sumw2 Analysis histos - simu and reco energy of triggered events
    h_triggered_BGOrec_E.Sumw2();
    h_triggered_BGOrec_E_corr.Sumw2();
    
    // Sumw2 Analysis histos - Geo
    h_geo_BGOrec_slopeX.Sumw2();
    h_geo_BGOrec_slopeY.Sumw2();
    h_geo_BGOrec_interceptX.Sumw2();
    h_geo_BGOrec_interceptY.Sumw2();
    h_geo_BGOreco_topMap.Sumw2();
    h_geo_BGOreco_bottomMap.Sumw2();
    h_layer_max_energy_ratio.Sumw2();

    // Sumw2 XTRL histos
    h_xtrl_energy_int.Sumw2();
    h_xtrl.Sumw2();

    // Create and load acceptance events cuts from config file
    cuts_conf flux_cuts;
    data_active_cuts active_cuts;
    load_flux_struct(flux_cuts, active_cuts, wd);

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
        
        // Fill the energy histos
        h_BGOrec_E.Fill(bgoTotalE_raw * _GeV);
        h_BGOrec_E_corr.Fill(bgoTotalE * _GeV);
        
        // Don't accept events outside the selected energy window
        if (bgoTotalE * _GeV < flux_cuts.min_event_energy || bgoTotalE * _GeV > flux_cuts.max_event_energy)
            continue;

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
            if (checkBGOreco_data(bgorec))
                h_trigger.Fill(bgoTotalE * _GeV);
            else
                continue;
        }
        else
            continue;
        
        best_track event_best_track;

        // Fill energy histos
        h_triggered_BGOrec_E.Fill(bgoTotalE_raw * _GeV);
        h_triggered_BGOrec_E_corr.Fill(bgoTotalE * _GeV);

        // Load BGO event class
        DmpBgoContainer bgoVault(DAMPE_bgo_nLayers);
        bgoVault.scanBGOHits(
            bgohits,
            bgoTotalE,
            DAMPE_bgo_nLayers);

        // evaluate the energy raio on each single layer of the BGO
        evaluateEnergyRatio(
            bgorec,
            flux_cuts,
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
        filter_geometric_cut = geometric_cut_data(bgorec);

        // **** BGO Fiducial Volume ****
        if (active_cuts.BGO_fiducial)
        {
            // maxElayer_cut
            filter_BGO_fiducial_maxElayer_cut = maxElayer_cut(
                bgorec,
                flux_cuts,
                bgoTotalE);

            // maxBarLayer_cut
            filter_BGO_fiducial_maxBarLayer_cut = maxBarLayer_cut(
                bgoVault.GetLayerBarNumber(),
                bgoVault.GetiMaxLayer(),
                bgoVault.GetIdxBarMaxLayer());

            // BGOTrackContainment_cut
            filter_BGO_fiducial_BGOTrackContainment_cut = BGOTrackContainment_cut(
                bgorec,
                flux_cuts);

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
                flux_cuts);
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
                    flux_cuts,
                    event_best_track);
            filter_all_cut *= filter_track_selection_cut;
        }

        // **** xtrl cut ****
        if (active_cuts.xtrl)
        {
            filter_xtrl_cut = xtrl_cut(
                bgoVault.GetSumRMS(),
                bgoVault.GetFracLayer(),
                flux_cuts,
                bgoTotalE,
                h_xtrl_energy_int,
                h_xtrl);
            filter_all_cut *= filter_xtrl_cut;
        }

        // **** psd_charge cut ****
        if (active_cuts.psd_charge)
        {
            filter_psd_charge_cut = psd_charge_cut(
                psdhits,
                bgorec,
                flux_cuts,
                event_best_track);
            filter_all_cut *= filter_psd_charge_cut;
        }

        // Fill cuts histos

        // Fill geometric cut histos
        if (filter_geometric_cut)
        {
            h_gometric_cut.Fill(bgoTotalE * _GeV);

            // Evaluate the position on the First BGO layer (after geometric cut)
            evaluateTopBottomPosition_data(
                bgorec,
                h_geo_BGOrec_slopeX,
                h_geo_BGOrec_slopeY,
                h_geo_BGOrec_interceptX,
                h_geo_BGOrec_interceptY,
                h_geo_BGOreco_topMap,
                h_geo_BGOreco_bottomMap);

            // Geometric cut && maxElayer cut
            if (filter_BGO_fiducial_maxElayer_cut)
                h_geometric_maxElayer_cut.Fill(bgoTotalE * _GeV);

            // Geometric cut && maxBarLayer cut
            if (filter_BGO_fiducial_maxBarLayer_cut)
                h_geometric_maxBarLayer_cut.Fill(bgoTotalE * _GeV);
            
            // Geometric cut && BGOTrackContainment cut
            if (filter_BGO_fiducial_BGOTrackContainment_cut)
                h_geometric_BGOTrackContainment_cut.Fill(bgoTotalE * _GeV);

            // Geometric cut && BGO fiducial cut
            if (filter_BGO_fiducial_cut)
                h_geometric_BGO_fiducial.Fill(bgoTotalE * _GeV);

            // Geometric cut && nBarLayer13 cut
            if (filter_nBarLayer13_cut)
                h_geometric_nBarLayer13_cut.Fill(bgoTotalE * _GeV);

            // Geometric cut && maxRms cut
            if (filter_maxRms_cut)
                h_geometric_maxRms_cut.Fill(bgoTotalE * _GeV);

            // Geometric cut && track selection cut
            if (filter_track_selection_cut)
                h_geometric_track_selection_cut.Fill(bgoTotalE * _GeV);

            // Geometric cut && XTRL cut
            if (filter_xtrl_cut)
                h_geometric_xtrl_cut.Fill(bgoTotalE * _GeV);

            // Geometric cut && PSD charge cut
            if (filter_psd_charge_cut)
                h_geometric_psd_charge_cut.Fill(bgoTotalE * _GeV);

            // Geometric cut and all cuts
            if (filter_all_cut)
                h_geometric_all_cut.Fill(bgoTotalE * _GeV);
        }

        // Fill BGO_fiducial_maxElayer cut histo
        if (filter_BGO_fiducial_maxElayer_cut)
            h_maxElayer_cut.Fill(bgoTotalE * _GeV);

        // Fill BGO_fiducial_maxBarLayer cut histo
        if (filter_BGO_fiducial_maxBarLayer_cut)
            h_maxBarLayer_cut.Fill(bgoTotalE * _GeV);

        // Fill BGO_fiducial_BGOTrackContainment cut histo
        if (filter_BGO_fiducial_BGOTrackContainment_cut)
            h_BGOTrackContainment_cut.Fill(bgoTotalE * _GeV);

        // Fill BGO fiducial volume cut
        if (filter_BGO_fiducial_cut)
        {
            h_BGO_fiducial_cut.Fill(bgoTotalE * _GeV);
        
            // BGO fiducial cut && nBarLayer13 cut
            if (filter_nBarLayer13_cut)
                h_BGOfiducial_nBarLayer13_cut.Fill(bgoTotalE * _GeV);
            
            // BGO fiducial cut && maxRms cut
            if (filter_maxRms_cut)
                h_BGOfiducial_maxRms_cut.Fill(bgoTotalE * _GeV);

            // BGO fiducial cut && track selection cut
            if (filter_track_selection_cut)
                h_BGOfiducial_track_selection_cut.Fill(bgoTotalE * _GeV);

            // BGO fiducial cut && XTRL cut
            if (filter_xtrl_cut)
                h_BGOfiducial_xtrl_cut.Fill(bgoTotalE * _GeV);

            // BGO fiducial cut && PSD charge cut
            if (filter_psd_charge_cut)
                h_BGOfiducial_psd_charge_cut.Fill(bgoTotalE * _GeV);

            // BGO fiducial cut && all cut
            if (filter_all_cut)
                h_BGOfiducial_all_cut.Fill(bgoTotalE * _GeV);
        }

        // Fill nBarLayer13 cut histo
        if (filter_nBarLayer13_cut)
            h_nBarLayer13_cut.Fill(bgoTotalE * _GeV);

        // Fill maxRms cut histo
        if (filter_maxRms_cut)
            h_maxRms_cut.Fill(bgoTotalE * _GeV);

        // Fill track selection cut histo
        if (filter_track_selection_cut)
            h_track_selection_cut.Fill(bgoTotalE * _GeV);

        // Fill XTRL cut histo
        if (filter_xtrl_cut)
            h_xtrl_cut.Fill(bgoTotalE * _GeV);

        // Fill PSD charge cut histo
        if (filter_psd_charge_cut)
            h_psd_charge_cut.Fill(bgoTotalE * _GeV);

        // Fill all cut histo
        if (active_cuts.nActiveCuts)
            if (filter_all_cut)
                h_all_cut.Fill(bgoTotalE * _GeV);
    }

    if (verbose)
    {
        std::cout << "\n\n ****** \n\n";
        std::cout << "Triggered events: " << h_trigger.GetEntries() << std::endl;

        auto refEntries = h_trigger.GetEntries();
        
        if (h_gometric_cut.GetEntries())
            std::cout << "geometric filtered events: " << h_gometric_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_gometric_cut.GetEntries()) / refEntries << std::endl;

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

    // Write histos to file
    // Acceptance - First-Cut histos
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

    // Create output ratio dir in the output TFile
    auto ratioDir = outFile.mkdir("Efficiency");

    // Create trigger folder
    auto trigger_dir = ratioDir->mkdir("Trigger");
    trigger_dir->cd();

    // Building ratio histos
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

    auto geo_analysisDir = outFile.mkdir("Analysis_GeoCut");
    geo_analysisDir->cd();
    
    h_geo_BGOrec_slopeX.Write();
    h_geo_BGOrec_slopeY.Write();
    h_geo_BGOrec_interceptX.Write();
    h_geo_BGOrec_interceptY.Write();
    h_geo_BGOreco_topMap.Write();
    h_geo_BGOreco_bottomMap.Write();

    outFile.cd();

    auto BGOdir = outFile.mkdir("BGO_Energy");
    BGOdir->cd();

    h_BGOrec_E.Write();
    h_BGOrec_E_corr.Write();
    h_layer_max_energy_ratio.Write();

    h_triggered_BGOrec_E.Write();
    h_triggered_BGOrec_E_corr.Write();
    
    for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
        h_layer_energy_ratio[lIdx].Write();

    auto XTRLdir = outFile.mkdir("xtrl");
    XTRLdir->cd();

    h_xtrl_energy_int.Write();
    h_xtrl.Write();

    outFile.cd();

    // Return all-cut histos
    return h_all_cut;
}