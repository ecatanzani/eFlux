#include "myHeader.h"
#include "acceptance.h"
#include "acceptance_cuts.h"
#include "energyMatch.h"

#include "TClonesArray.h"
#include "TString.h"

/**
 * @brief 
 * 
 * @param accInputPath 
 * @param verbose 
 */

std::string getListPath(const std::string accInputPath, const bool MC)
{
    std::string MClist = accInputPath;
    std::string relMClist;
    if (MC)
    {
        std::string relMClist = "/MC.txt";
        MClist.append(relMClist);
    }
    else
    {
        std::string relMClist = "/Data.txt";
        MClist.append(relMClist);
    }
    return MClist;
}

std::shared_ptr<DmpChain> aggregateEventsDmpChain(
    const std::string accInputPath,
    const bool verbose)
{
    // ****** Access data using DAMPE Chain ******

    // Create DmpChain object
    std::shared_ptr<DmpChain> dmpch = std::make_shared<DmpChain>("CollectionTree");

    // Add MC file list to DmpChain
    //dmpch->AddFromList(getListPath(accInputPath, true).c_str());
    dmpch->AddFromList(accInputPath.c_str());
    if (verbose)
        dmpch->GetListOfFiles()->Print();

    return dmpch;
}

std::shared_ptr<TChain> aggregateEventsTChain(
    const std::string accInputPath,
    const bool verbose)
{
    // ****** Access data using ROOT TChain ******

    // Create TChain object
    //TChain* dmpch = new TChain("CollectionTree");
    std::shared_ptr<TChain> dmpch = std::make_shared<TChain>("CollectionTree");
    //std::shared_ptr<TChain> dmpch( new TChain("CollectionTree") );

    // Reading list of MC files
    //std::ifstream input_file(getListPath(accInputPath, true).c_str());
    std::ifstream input_file(accInputPath.c_str());
    if (!input_file.is_open())
    {
        std::cerr << "\nERROR 100! File not open " << getListPath(accInputPath, true) << "\n\n";
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    while (input_stream >> tmp_str)
    {
        dmpch->Add(tmp_str.c_str());
        if (verbose)
            std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }

    return dmpch;
}

double wtsydp(
    const float minene,
    const float maxene,
    const float index)
{
    float dene = maxene - minene;
    if (index != -1)
        return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * dene)), 1. / index);
    else
        return dene / log(maxene / minene);
}

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

    // **** Acceptance

    // First-Cut histos
    TH1D h_incoming("h_incoming", "Energy Distribution of the incoming particles", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_trigger("h_trigger", "Energy Distribution of the triggered particles", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_gometric_cut("h_gometric_cut", "Energy Distribution - geometric cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxElayer_cut("h_maxElayer_cut", "Energy Distribution - maxElayer cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxBarLayer_cut("h_maxBarLayer_cut", "Energy Distribution - maxBarLayer cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOTrackContainment_cut("h_BGOTrackContainment_cut", "Energy Distribution - BGOTrackContainment cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGO_fiducial("h_BGO_fiducial", "Energy Distibution - BGO fiducial cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_nBarLayer13_cut("h_nBarLayer13_cut", "Energy Distribution - nBarLayer13 cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxRms_cut("h_maxRms_cut", "Energy Distribution - maxRms cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_track_selection_cut("h_track_selection_cut", "Energy Distribution - track selection cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_xtrl_cut("h_xtrl_cut", "Energy Distribution - xtrl cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_psd_charge_cut("h_psd_charge_cut", "Energy Distribution - psd charge cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_all_cut("h_all_cut", "Energy Distribution - All cut ", logEBins.size() - 1, &(logEBins[0]));

    // **** Analysis histos

    TH1D h_BGOrec_E("h_BGOrec_E", "BGO Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOrec_E_corr("h_BGOrec_E_corr", "BGO Corrected Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_simu_energy("h_simu_energy", "Simu Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_energy_diff("h_energy_diff", "Simu vs Corrected Reco BGO energy", 50, -100, 100);

    TH1D h_accepted_BGOrec_E("h_accepted_BGOrec_E", "BGO Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_accepted_BGOrec_E_corr("h_accepted_BGOrec_E_corr", "BGO Corrected Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_accepted_simu_energy("h_accepted_simu_energy", "Simu Energy", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_accepted_energy_diff("h_accepted_energy_diff", "Simu vs Corrected Reco BGO energy", 50, -100, 100);

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

    // Ratio of layer energy respect to total BGO energy
    TH1D h_preGeo_layer_energy_ratio("h_preGeo_layer_energy_ratio", "Layer Energy Ratio", 100, 0, 1);

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

    // *****
    h_BGOrec_E.Sumw2();
    h_BGOrec_E_corr.Sumw2();
    h_simu_energy.Sumw2();
    h_energy_diff.Sumw2();

    h_accepted_BGOrec_E.Sumw2();
    h_accepted_BGOrec_E_corr.Sumw2();
    h_accepted_simu_energy.Sumw2();
    h_accepted_energy_diff.Sumw2();

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

    h_incoming.Sumw2();
    h_trigger.Sumw2();
    h_gometric_cut.Sumw2();
    h_maxElayer_cut.Sumw2();
    h_maxBarLayer_cut.Sumw2();
    h_BGOTrackContainment_cut.Sumw2();
    h_BGO_fiducial.Sumw2();
    h_nBarLayer13_cut.Sumw2();
    h_maxRms_cut.Sumw2();
    h_track_selection_cut.Sumw2();
    h_xtrl_cut.Sumw2();
    h_psd_charge_cut.Sumw2();
    h_all_cut.Sumw2();

    h_preGeo_layer_energy_ratio.Sumw2();
    h_layer_max_energy_ratio.Sumw2();

    // Create and load acceptance events cuts from config file
    acceptance_conf acceptance_cuts;
    acceptance_active_cuts active_cuts;
    load_acceptance_struct(acceptance_cuts, active_cuts, wd);

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
            h_trigger.Fill(simuEnergy * _GeV);
            if (checkBGOreco(bgorec, simu_primaries))
                h_incoming.Fill(simuEnergy * _GeV);
            else
                continue;
        }
        else
        {
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

        // Fill geometric-cut histo
        if (active_cuts.geometry)
            if (geometric_cut(simu_primaries))
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
            }
        
        bool filter_nBarLayer13_cut = false;
        bool filter_maxRms_cut = false;
        bool filter_track_selection_cut = false;
        bool filter_xtrl_cut = false;
        bool filter_psd_charge_cut = false;

        best_track event_best_track;
        bool all_event_filter = true;

        // Fill energy histos
        h_accepted_BGOrec_E.Fill(bgoTotalE_raw * _GeV);
        h_accepted_BGOrec_E_corr.Fill(bgoTotalE * _GeV);
        h_accepted_simu_energy.Fill(simuEnergy * _GeV);
        h_accepted_energy_diff.Fill((simuEnergy - bgoTotalE) * _GeV);

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

        // maxElayer_cut && geometric cut
        auto filter_maxElayer_cut = maxElayer_cut(
            bgorec,
            acceptance_cuts,
            bgoTotalE);

        // maxBarLayer_cut && geometric cut
        auto filter_maxBarLayer_cut = maxBarLayer_cut(
            bgoVault.GetLayerBarNumber(),
            bgoVault.GetiMaxLayer(),
            bgoVault.GetIdxBarMaxLayer());

        // BGOTrackContainment_cut && geometric cut
        auto filter_BGOTrackContainment_cut = BGOTrackContainment_cut(
            bgorec,
            acceptance_cuts);

        if (filter_maxElayer_cut)
            h_maxElayer_cut.Fill(simuEnergy * _GeV);
        if (filter_maxBarLayer_cut)
            h_maxBarLayer_cut.Fill(simuEnergy * _GeV);
        if (filter_BGOTrackContainment_cut)
            h_BGOTrackContainment_cut.Fill(simuEnergy * _GeV);

        // **** BGO Fiducial Volume ****
        if (active_cuts.BGO_fiducial)
        {
            bool filter_BGO_fiducial = true;

            filter_BGO_fiducial *= filter_maxElayer_cut;
            filter_BGO_fiducial *= filter_maxBarLayer_cut;
            filter_BGO_fiducial *= filter_BGOTrackContainment_cut;

            // BGO_fiducial_cut
            if (filter_BGO_fiducial)
                h_BGO_fiducial.Fill(simuEnergy * _GeV);
        }

        // **** nBarLayer13 cut ****
        if (active_cuts.nBarLayer13)
        {
            filter_nBarLayer13_cut = nBarLayer13_cut(
                bgohits,
                bgoVault.GetSingleLayerBarNumber(13),
                bgoTotalE);
            all_event_filter *= filter_nBarLayer13_cut;
            if (filter_nBarLayer13_cut)
                h_nBarLayer13_cut.Fill(simuEnergy * _GeV);
        }

        // **** maxRms cut ****
        if (active_cuts.maxRms)
        {
            filter_maxRms_cut = maxRms_cut(
                bgoVault.GetLayerBarNumber(),
                bgoVault.GetRmsLayer(),
                bgoTotalE,
                acceptance_cuts);
            all_event_filter *= filter_maxRms_cut;
            if (filter_maxRms_cut)
                h_maxRms_cut.Fill(simuEnergy * _GeV);
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
            all_event_filter *= filter_track_selection_cut;
            if (filter_track_selection_cut)
                h_track_selection_cut.Fill(simuEnergy * _GeV);
        }

        // **** xtrl cut ****
        if (active_cuts.xtrl)
        {
            filter_xtrl_cut = xtrl_cut(
                bgoVault.GetSumRMS(),
                bgoVault.GetFracLayer(),
                acceptance_cuts);
            all_event_filter *= filter_xtrl_cut;
            if (filter_xtrl_cut)
                h_xtrl_cut.Fill(simuEnergy * _GeV);
        }

        // **** psd_charge cut ****
        if (active_cuts.psd_charge)
        {
            filter_psd_charge_cut = psd_charge_cut(
                psdhits,
                bgorec,
                acceptance_cuts,
                event_best_track);
            all_event_filter *= filter_psd_charge_cut;
            if (filter_psd_charge_cut)
                h_psd_charge_cut.Fill(simuEnergy * _GeV);
        }

        // **** All-cuts ****
        if (active_cuts.nActiveCuts)
            if (all_event_filter)
                h_all_cut.Fill(simuEnergy * _GeV);
    }

    if (verbose)
    {
        std::cout << "\n\n ****** \n\n";
        std::cout << "generated events in good energy range: " << h_incoming.GetEntries() << std::endl;
        std::cout << "triggered events: " << h_trigger.GetEntries() << std::endl;

        auto refEntries = h_trigger.GetEntries();

        if (h_gometric_cut.GetEntries())
            std::cout << "geometric filtered events: " << h_gometric_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_gometric_cut.GetEntries()) / refEntries << std::endl;

        if (h_BGO_fiducial.GetEntries())
            std::cout << "BGO fiducial filtered events: " << h_BGO_fiducial.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_BGO_fiducial.GetEntries()) / refEntries << std::endl;

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
            std::cout << "psd_charge filtered events: " << h_all_cut.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_all_cut.GetEntries()) / refEntries;

        std::cout << "\n\n ****** \n\n";
    }

    //double genSurface = 4 * TMath::Pi() * pow(acceptance_cuts.vertex_radius, 2) / 2;
    double genSurface = 8 * pow(TMath::Pi(),2) * pow(acceptance_cuts.vertex_radius, 2) / 2;

    // Building acceptance histos
    auto h_acceptance_gometric_cut = static_cast<TH1D *>(h_gometric_cut.Clone("h_acceptance_gometric_cut"));
    auto h_acceptance_maxElayer_cut = static_cast<TH1D *>(h_maxElayer_cut.Clone("h_acceptance_maxElayer_cut"));
    auto h_acceptance_maxBarLayer_cut = static_cast<TH1D *>(h_maxBarLayer_cut.Clone("h_acceptance_maxBarLayer_cut"));
    auto h_acceptance_BGOTrackContainment_cut = static_cast<TH1D *>(h_BGOTrackContainment_cut.Clone("h_acceptance_BGOTrackContainment_cut"));
    auto h_acceptance_BGO_fiducial = static_cast<TH1D *>(h_BGO_fiducial.Clone("h_acceptance_BGO_fiducial"));
    auto h_acceptance_nBarLayer13_cut = static_cast<TH1D *>(h_nBarLayer13_cut.Clone("h_acceptance_nBarLayer13_cut"));
    auto h_acceptance_maxRms_cut = static_cast<TH1D *>(h_maxRms_cut.Clone("h_acceptance_maxRms_cut"));
    auto h_acceptance_track_selection_cut = static_cast<TH1D *>(h_track_selection_cut.Clone("h_acceptance_track_selection_cut"));
    auto h_acceptance_xtrl_cut = static_cast<TH1D *>(h_xtrl_cut.Clone("h_acceptance_xtrl_cut"));
    auto h_acceptance_psd_charge_cut = static_cast<TH1D *>(h_psd_charge_cut.Clone("h_acceptance_psd_charge_cut"));
    auto h_acceptance_all_cut = static_cast<TH1D *>(h_all_cut.Clone("h_acceptance_all_cut"));

    h_acceptance_gometric_cut->Divide(&h_incoming);
    h_acceptance_maxElayer_cut->Divide(&h_incoming);
    h_acceptance_maxBarLayer_cut->Divide(&h_incoming);
    h_acceptance_BGOTrackContainment_cut->Divide(&h_incoming);
    h_acceptance_BGO_fiducial->Divide(&h_incoming);
    h_acceptance_nBarLayer13_cut->Divide(&h_incoming);
    h_acceptance_maxRms_cut->Divide(&h_incoming);
    h_acceptance_track_selection_cut->Divide(&h_incoming);
    h_acceptance_xtrl_cut->Divide(&h_incoming);
    h_acceptance_psd_charge_cut->Divide(&h_incoming);
    h_acceptance_all_cut->Divide(&h_incoming);

    h_acceptance_gometric_cut->Scale(genSurface);
    h_acceptance_maxElayer_cut->Scale(genSurface);
    h_acceptance_maxBarLayer_cut->Scale(genSurface);
    h_acceptance_BGOTrackContainment_cut->Scale(genSurface);
    h_acceptance_BGO_fiducial->Scale(genSurface);
    h_acceptance_nBarLayer13_cut->Scale(genSurface);
    h_acceptance_maxRms_cut->Scale(genSurface);
    h_acceptance_track_selection_cut->Scale(genSurface);
    h_acceptance_xtrl_cut->Scale(genSurface);
    h_acceptance_psd_charge_cut->Scale(genSurface);
    h_acceptance_all_cut->Scale(genSurface);

    // Builing vectors
    std::vector<double> energyValues(h_incoming.GetXaxis()->GetNbins(), 0);

    std::vector<double> acceptanceValues_gometric_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxElayer_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxBarLayer_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGOTrackContainment_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGO_fiducial_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_nBarLayer13_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxRms_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_track_selection_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_xtrl_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_psd_charge_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_all_cut(energyValues.size(), 0);

    for (auto it = logEBins.begin(); it != (logEBins.end() - 1); ++it)
    {
        auto index = std::distance(logEBins.begin(), it);
        energyValues[index] = wtsydp(*it, *(it + 1), -1);
        acceptanceValues_gometric_cut[index] = h_acceptance_gometric_cut->GetBinContent(index + 1);
        acceptanceValues_maxElayer_cut[index] = h_acceptance_maxElayer_cut->GetBinContent(index + 1);
        acceptanceValues_maxBarLayer_cut[index] = h_acceptance_maxBarLayer_cut->GetBinContent(index + 1);
        acceptanceValues_BGOTrackContainment_cut[index] = h_acceptance_BGOTrackContainment_cut->GetBinContent(index + 1);
        acceptanceValues_BGO_fiducial_cut[index] = h_acceptance_BGO_fiducial->GetBinContent(index + 1);
        acceptanceValues_nBarLayer13_cut[index] = h_acceptance_nBarLayer13_cut->GetBinContent(index + 1);
        acceptanceValues_maxRms_cut[index] = h_acceptance_maxRms_cut->GetBinContent(index + 1);
        acceptanceValues_track_selection_cut[index] = h_acceptance_track_selection_cut->GetBinContent(index + 1);
        acceptanceValues_xtrl_cut[index] = h_acceptance_xtrl_cut->GetBinContent(index + 1);
        acceptanceValues_psd_charge_cut[index] = h_acceptance_psd_charge_cut->GetBinContent(index + 1);
        acceptanceValues_all_cut[index] = h_acceptance_all_cut->GetBinContent(index + 1);
    }

    // Building graphs
    TGraph gr_acceptance_gometric_cut(energyValues.size(), &energyValues[0], &acceptanceValues_gometric_cut[0]);
    TGraph gr_acceptance_maxElayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxElayer_cut[0]);
    TGraph gr_acceptance_maxBarLayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxBarLayer_cut[0]);
    TGraph gr_acceptance_BGOTrackContainment_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGOTrackContainment_cut[0]);
    TGraph gr_acceptance_BGO_fiducial_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_cut[0]);
    TGraph gr_acceptance_nBarLayer13_cut(energyValues.size(), &energyValues[0], &acceptanceValues_nBarLayer13_cut[0]);
    TGraph gr_acceptance_maxRms_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxRms_cut[0]);
    TGraph gr_acceptance_track_selection_cut(energyValues.size(), &energyValues[0], &acceptanceValues_track_selection_cut[0]);
    TGraph gr_acceptance_xtrl_cut(energyValues.size(), &energyValues[0], &acceptanceValues_xtrl_cut[0]);
    TGraph gr_acceptance_psd_charge_cut(energyValues.size(), &energyValues[0], &acceptanceValues_psd_charge_cut[0]);
    TGraph gr_acceptance_all_cut(energyValues.size(), &energyValues[0], &acceptanceValues_all_cut[0]);

    gr_acceptance_gometric_cut.SetName("gr_acceptance_gometric_cut");
    gr_acceptance_maxElayer_cut.SetName("gr_acceptance_maxElayer_cut");
    gr_acceptance_maxBarLayer_cut.SetName("gr_acceptance_maxBarLayer_cut");
    gr_acceptance_BGOTrackContainment_cut.SetName("gr_acceptance_BGOTrackContainment_cut");
    gr_acceptance_BGO_fiducial_cut.SetName("gr_acceptance_BGO_fiducial_cut");
    gr_acceptance_nBarLayer13_cut.SetName("gr_acceptance_nBarLayer13_cut");
    gr_acceptance_maxRms_cut.SetName("gr_acceptance_maxRms_cut");
    gr_acceptance_track_selection_cut.SetName("gr_acceptance_track_selection_cut");
    gr_acceptance_xtrl_cut.SetName("gr_acceptance_xtrl_cut");
    gr_acceptance_psd_charge_cut.SetName("gr_acceptance_psd_charge_cut");
    gr_acceptance_all_cut.SetName("gr_acceptance_all_cut");

    gr_acceptance_gometric_cut.SetTitle("Acceptance - geometric cut");
    gr_acceptance_maxElayer_cut.SetTitle("Acceptance - maxElateral cut");
    gr_acceptance_maxBarLayer_cut.SetTitle("Acceptance - maxBarLayer cut");
    gr_acceptance_BGOTrackContainment_cut.SetTitle("Acceptance - BGOTrackContainment cut");
    gr_acceptance_BGO_fiducial_cut.SetTitle("Acceptance - BGO fiducial volume cut");
    gr_acceptance_nBarLayer13_cut.SetTitle("Acceptance - nBarLayer13 cut");
    gr_acceptance_maxRms_cut.SetTitle("Acceptance - maxRms cut");
    gr_acceptance_track_selection_cut.SetTitle("Acceptance - track selection cut");
    gr_acceptance_xtrl_cut.SetTitle("Acceptance - XTRL cut");
    gr_acceptance_psd_charge_cut.SetTitle("Acceptance - PSD charge selection cut");
    gr_acceptance_all_cut.SetTitle("Acceptance - all cut");

    // Write histos to file
    h_incoming.Write();
    h_trigger.Write();
    h_gometric_cut.Write();
    h_maxElayer_cut.Write();
    h_maxBarLayer_cut.Write();
    h_BGOTrackContainment_cut.Write();
    h_BGO_fiducial.Write();
    h_nBarLayer13_cut.Write();
    h_maxRms_cut.Write();
    h_track_selection_cut.Write();
    h_xtrl_cut.Write();
    h_psd_charge_cut.Write();
    h_all_cut.Write();

    // Create output acceptance dir in the output TFile
    auto acceptanceDir = outFile.mkdir("Acceptance");
    acceptanceDir->cd();

    // Write final TGraphs
    gr_acceptance_gometric_cut.Write();
    gr_acceptance_maxElayer_cut.Write();
    gr_acceptance_maxBarLayer_cut.Write();
    gr_acceptance_BGOTrackContainment_cut.Write();
    gr_acceptance_BGO_fiducial_cut.Write();
    gr_acceptance_nBarLayer13_cut.Write();
    gr_acceptance_maxRms_cut.Write();
    gr_acceptance_track_selection_cut.Write();
    gr_acceptance_xtrl_cut.Write();
    gr_acceptance_psd_charge_cut.Write();
    gr_acceptance_all_cut.Write();

    // Return to main TFile directory
    outFile.cd();

    // Create output ratio dir in the output TFile
    auto ratioDir = outFile.mkdir("Ratios");
    ratioDir->cd();

    // Building ratio histos
    auto h_ratio_gometric_cut = static_cast<TH1D *>(h_gometric_cut.Clone("h_ratio_gometric_cut"));
    auto h_ratio_maxElayer_cut = static_cast<TH1D *>(h_maxElayer_cut.Clone("h_ratio_maxElayer_cut"));
    auto h_ratio_maxBarLayer_cut = static_cast<TH1D *>(h_maxBarLayer_cut.Clone("h_ratio_maxBarLayer_cut"));
    auto h_ratio_BGOTrackContainment_cut = static_cast<TH1D *>(h_BGOTrackContainment_cut.Clone("h_ratio_BGOTrackContainment_cut"));
    auto h_ratio_BGO_fiducial = static_cast<TH1D *>(h_BGO_fiducial.Clone("h_ratio_BGO_fiducial"));
    auto h_ratio_nBarLayer13_cut = static_cast<TH1D *>(h_nBarLayer13_cut.Clone("h_ratio_nBarLayer13_cut"));
    auto h_ratio_maxRms_cut = static_cast<TH1D *>(h_maxRms_cut.Clone("h_ratio_maxRms_cut"));
    auto h_ratio_track_selection_cut = static_cast<TH1D *>(h_track_selection_cut.Clone("h_ratio_track_selection_cut"));
    auto h_ratio_xtrl_cut = static_cast<TH1D *>(h_xtrl_cut.Clone("h_ratio_xtrl_cut"));
    auto h_ratio_psd_charge_cut = static_cast<TH1D *>(h_psd_charge_cut.Clone("h_ratio_psd_charge_cut"));
    auto h_ratio_all_cut = static_cast<TH1D *>(h_all_cut.Clone("h_ratio_all_cut"));

    h_ratio_gometric_cut->Divide(&h_trigger);
    h_ratio_maxElayer_cut->Divide(&h_trigger);
    h_ratio_maxBarLayer_cut->Divide(&h_trigger);
    h_ratio_BGOTrackContainment_cut->Divide(&h_trigger);
    h_ratio_BGO_fiducial->Divide(&h_trigger);
    h_ratio_nBarLayer13_cut->Divide(&h_trigger);
    h_ratio_maxRms_cut->Divide(&h_trigger);
    h_ratio_track_selection_cut->Divide(&h_trigger);
    h_ratio_xtrl_cut->Divide(&h_trigger);
    h_ratio_psd_charge_cut->Divide(&h_trigger);
    h_ratio_all_cut->Divide(&h_trigger);

    h_ratio_gometric_cut->Write();
    h_ratio_maxElayer_cut->Write();
    h_ratio_maxBarLayer_cut->Write();
    h_ratio_BGOTrackContainment_cut->Write();
    h_ratio_BGO_fiducial->Write();
    h_ratio_nBarLayer13_cut->Write();
    h_ratio_maxRms_cut->Write();
    h_ratio_track_selection_cut->Write();
    h_ratio_xtrl_cut->Write();
    h_ratio_psd_charge_cut->Write();
    h_ratio_all_cut->Write();

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

    h_preGeo_layer_energy_ratio.Write();

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

    h_accepted_BGOrec_E.Write();
    h_accepted_BGOrec_E_corr.Write();
    h_accepted_simu_energy.Write();
    h_accepted_energy_diff.Write();

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