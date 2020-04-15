/*
    Consider this repo as a manual: https://github.com/DAMPEEU/DmpTools/blob/master/HE_skimmer/code/app/main.cc#L318
*/

#include "myHeader.h"
#include "acceptance.h"
#include "acceptance_cuts.h"
#include "energyMatch.h"

#include "TClonesArray.h"

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

inline std::shared_ptr<TH1D> buildHistoFromVector(
    const std::vector<double> &energyValues,
    const std::vector<double> &consgFactor)
{
    std::shared_ptr<TH1D> histo = std::make_shared<TH1D>("histo", "histoTitle", consgFactor.size(), energyValues[0], energyValues[energyValues.size() - 1]);
    for (auto idx = 1; idx <= histo->GetNbinsX(); ++idx)
        histo->SetBinContent(idx, consgFactor[idx - 1]);

    return histo;
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

    // Cut histos
    TH1D h_incoming("h_incoming", "Energy Distribution of the incoming particles", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_gometric_cut("h_gometric_cut", "Energy Distribution - geometric cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxElayer_cut("h_maxElayer_cut", "Energy Distribution - maxElateral cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxBarLayer_cut("h_maxBarLayer_cut", "Energy Distribution - maxBarLayer cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_BGOTrackContainment_cut("h_BGOTrackContainment_cut", "Energy Distribution - BGOTrackContainment cut ", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_nBarLayer13_cut("h_nBarLayer13_cut", "Energy Distribution - nBarLayer13 cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_maxRms_cut("h_maxRms_cut", "Energy Distribution - maxRms cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_track_selection_cut("h_track_selection_cut", "Energy Distribution - track selection cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_xtrl_cut("h_xtrl_cut", "Energy Distribution - xtrl cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_psd_charge_cut("h_psd_charge_cut", "Energy Distribution - psd charge cut", logEBins.size() - 1, &(logEBins[0]));
    TH1D h_all_cut("h_all_cut", "Energy Distribution - All cut ", logEBins.size() - 1, &(logEBins[0]));

    h_incoming.Sumw2();
    h_gometric_cut.Sumw2();
    h_maxElayer_cut.Sumw2();
    h_maxBarLayer_cut.Sumw2();
    h_BGOTrackContainment_cut.Sumw2();
    h_nBarLayer13_cut.Sumw2();
    h_maxRms_cut.Sumw2();
    h_track_selection_cut.Sumw2();
    h_xtrl_cut.Sumw2();
    h_psd_charge_cut.Sumw2();
    h_all_cut.Sumw2();

    // Create and load acceptance events cuts from config file
    acceptance_conf acceptance_cuts;
    acceptance_active_cuts active_cuts;
    load_acceptance_struct(acceptance_cuts, active_cuts, wd);

    double _GeV = 0.001;

    for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx)
    {
        // Get chain event
        dmpch->GetEvent(evIdx);

        // GOOD EVENT variable
        bool passEvent = true;

        // Event printout
        if (verbose)
            if (((evIdx + 1) % _kStep) == 0)
                std::cout << "\nProcessed " << evIdx + 1 << " events / " << nevents;

        // Get event total energy
        //double bgoTotalE = bgorec->GetTotalEnergy();      // Energy in MeV - not corrected
        double bgoTotalE = bgorec->GetElectronEcor();    // Returns corrected energy assuming this was an electron (MeV)
        double simuEnergy = simu_primaries->pvpart_ekin; //Energy of simu primaries particle in MeV

        // Don't accept events outside the selected energy window
        if (simuEnergy * _GeV < acceptance_cuts.min_event_energy || simuEnergy * _GeV > acceptance_cuts.max_event_energy)
            continue;

        h_incoming.Fill(simuEnergy * _GeV);

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

        bool filter_maxElayer_cut = false;
        bool filter_maxBarLayer_cut = false;
        bool filter_BGOTrackContainment_cut = false;
        bool filter_nBarLayer13_cut = false;
        bool filter_maxRms_cut = false;
        bool filter_track_selection_cut = false;
        bool filter_xtrl_cut = false;
        bool filter_psd_charge_cut = false;

        best_track event_best_track;
        bool all_event_filter = false;

        // **** First cuts ****
        
        if (geometric_cut(simu_primaries))
        {
            h_gometric_cut.Fill(simuEnergy * _GeV);
            if (active_cuts.maxElater)
            {
                filter_maxElayer_cut = maxElayer_cut(bgorec, acceptance_cuts, bgoTotalE);
                all_event_filter = filter_maxElayer_cut;
                if (filter_maxElayer_cut)
                    h_maxElayer_cut.Fill(simuEnergy * _GeV);
            }
            if (active_cuts.maxBarLayer)
            {
                filter_maxBarLayer_cut = maxBarLayer_cut(bgohits, nBgoHits);
                all_event_filter *= filter_maxBarLayer_cut;
                if (filter_maxBarLayer_cut)
                    h_maxBarLayer_cut.Fill(simuEnergy * _GeV);
            }
            if (active_cuts.BGOTrackContainment)
            {
                filter_BGOTrackContainment_cut = BGOTrackContainment_cut(bgorec, acceptance_cuts, passEvent);
                all_event_filter *= filter_BGOTrackContainment_cut;
                if (filter_BGOTrackContainment_cut)
                    h_BGOTrackContainment_cut.Fill(simuEnergy * _GeV);
            }
            if (active_cuts.nBarLayer13)
            {
                filter_nBarLayer13_cut = nBarLayer13_cut(bgohits, layerBarNumber[13], bgoTotalE);
                all_event_filter *= filter_nBarLayer13_cut;
                if (filter_nBarLayer13_cut)
                    h_nBarLayer13_cut.Fill(simuEnergy * _GeV);
            }
            if (active_cuts.maxRms)
            {
                filter_maxRms_cut = maxRms_cut(layerBarNumber, rmsLayer, bgoTotalE, acceptance_cuts);
                all_event_filter *= filter_maxRms_cut;
                if (filter_maxRms_cut)
                    h_maxRms_cut.Fill(simuEnergy * _GeV);
            }
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
            if (active_cuts.xtrl)
            {
                filter_xtrl_cut = xtrl_cut(sumRms, fracLayer, acceptance_cuts);
                all_event_filter *= filter_xtrl_cut;
                if (filter_xtrl_cut)
                    h_xtrl_cut.Fill(simuEnergy * _GeV);
            }
            if (active_cuts.psd_charge)
            {
                filter_psd_charge_cut = psd_charge_cut(psdhits, bgorec, acceptance_cuts, event_best_track);
                all_event_filter *= filter_psd_charge_cut;
                if (filter_psd_charge_cut)
                    h_psd_charge_cut.Fill(simuEnergy * _GeV);
            }

            // **** All-cuts ****
            if (all_event_filter)
                h_all_cut.Fill(simuEnergy * _GeV);
        }
    }

    if (verbose)
    {
        std::cout << "\n\ngeometric filtered events: " << h_gometric_cut.GetEntries() << "/" << nevents << std::endl;
        std::cout << std::endl << "all-cut filtered events: " << h_all_cut.GetEntries() << "/" << nevents << "\n\n";
    }

    double genSurface = 4 * TMath::Pi() * pow(acceptance_cuts.vertex_radius, 2) / 2;

    // Building acceptance histos
    auto h_acceptance_gometric_cut = static_cast<TH1D *>(h_gometric_cut.Clone("h_acceptance_gometric_cut"));
    auto h_acceptance_maxElateral_cut = static_cast<TH1D *>(h_maxElayer_cut.Clone("h_acceptance_maxElateral_cut"));
    auto h_acceptance_maxBarLayer_cut = static_cast<TH1D *>(h_maxBarLayer_cut.Clone("h_acceptance_maxBarLayer_cut"));
    auto h_acceptance_BGOTrackContainment_cut = static_cast<TH1D *>(h_BGOTrackContainment_cut.Clone("h_acceptance_BGOTrackContainment_cut"));
    auto h_acceptance_nBarLayer13_cut = static_cast<TH1D *>(h_nBarLayer13_cut.Clone("h_acceptance_nBarLayer13_cut"));
    auto h_acceptance_maxRms_cut = static_cast<TH1D *>(h_maxRms_cut.Clone("h_acceptance_maxRms_cut"));
    auto h_acceptance_track_selection_cut = static_cast<TH1D *>(h_track_selection_cut.Clone("h_acceptance_track_selection_cut"));
    auto h_acceptance_xtrl_cut = static_cast<TH1D *>(h_xtrl_cut.Clone("h_acceptance_xtrl_cut"));
    auto h_acceptance_psd_charge_cut = static_cast<TH1D *>(h_psd_charge_cut.Clone("h_acceptance_psd_charge_cut"));
    auto h_acceptance_all_cut = static_cast<TH1D *>(h_all_cut.Clone("h_acceptance_all_cut"));

    h_acceptance_gometric_cut->Scale(genSurface);
    h_acceptance_maxElateral_cut->Scale(genSurface);
    h_acceptance_maxBarLayer_cut->Scale(genSurface);
    h_acceptance_BGOTrackContainment_cut->Scale(genSurface);
    h_acceptance_nBarLayer13_cut->Scale(genSurface);
    h_acceptance_maxRms_cut->Scale(genSurface);
    h_acceptance_track_selection_cut->Scale(genSurface);
    h_acceptance_xtrl_cut->Scale(genSurface);
    h_acceptance_psd_charge_cut->Scale(genSurface);
    h_acceptance_all_cut->Scale(genSurface);

    h_acceptance_gometric_cut->Divide(&h_incoming);
    h_acceptance_maxElateral_cut->Divide(&h_incoming);
    h_acceptance_maxBarLayer_cut->Divide(&h_incoming);
    h_acceptance_BGOTrackContainment_cut->Divide(&h_incoming);
    h_acceptance_nBarLayer13_cut->Divide(&h_incoming);
    h_acceptance_maxRms_cut->Divide(&h_incoming);
    h_acceptance_track_selection_cut->Divide(&h_incoming);
    h_acceptance_xtrl_cut->Divide(&h_incoming);
    h_acceptance_psd_charge_cut->Divide(&h_incoming);
    h_acceptance_all_cut->Divide(&h_incoming);

    // Builing vectors
    std::vector<double> energyValues(h_acceptance_maxElateral_cut->GetXaxis()->GetNbins(), 0);

    std::vector<double> acceptanceValues_gometric_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxElateral_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxBarLayer_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGOTrackContainment_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_nBarLayer13_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxRms_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_track_selection_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_xtrl_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_psd_charge_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_all_cut(energyValues.size(), 0);

    for (auto it = logEBins.begin(); it != (logEBins.end()-1); ++it)
    {
        auto index = std::distance(logEBins.begin(), it);
        energyValues[index] = wtsydp(*it, *(it + 1), -1);
        acceptanceValues_gometric_cut[index] = h_acceptance_gometric_cut->GetBinContent(index + 1);
        acceptanceValues_maxElateral_cut[index] = h_acceptance_maxElateral_cut->GetBinContent(index + 1);
        acceptanceValues_maxBarLayer_cut[index] = h_acceptance_maxBarLayer_cut->GetBinContent(index + 1);
        acceptanceValues_BGOTrackContainment_cut[index] = h_acceptance_BGOTrackContainment_cut->GetBinContent(index + 1);
        acceptanceValues_nBarLayer13_cut[index] = h_acceptance_nBarLayer13_cut->GetBinContent(index + 1);
        acceptanceValues_maxRms_cut[index] = h_acceptance_maxRms_cut->GetBinContent(index + 1);
        acceptanceValues_track_selection_cut[index] = h_acceptance_track_selection_cut->GetBinContent(index + 1);
        acceptanceValues_xtrl_cut[index] = h_acceptance_xtrl_cut->GetBinContent(index + 1);
        acceptanceValues_psd_charge_cut[index] = h_acceptance_psd_charge_cut->GetBinContent(index + 1);
        acceptanceValues_all_cut[index] = h_acceptance_all_cut->GetBinContent(index + 1);
    }

    // Building graphs
    TGraph gr_acceptance_gometric_cut(energyValues.size(), &energyValues[0], &acceptanceValues_gometric_cut[0]);
    TGraph gr_acceptance_maxElateral_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxElateral_cut[0]);
    TGraph gr_acceptance_maxBarLayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxBarLayer_cut[0]);
    TGraph gr_acceptance_BGOTrackContainment_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGOTrackContainment_cut[0]);
    TGraph gr_acceptance_nBarLayer13_cut(energyValues.size(), &energyValues[0], &acceptanceValues_nBarLayer13_cut[0]);
    TGraph gr_acceptance_maxRms_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxRms_cut[0]);
    TGraph gr_acceptance_track_selection_cut(energyValues.size(), &energyValues[0], &acceptanceValues_track_selection_cut[0]);
    TGraph gr_acceptance_xtrl_cut(energyValues.size(), &energyValues[0], &acceptanceValues_xtrl_cut[0]);
    TGraph gr_acceptance_psd_charge_cut(energyValues.size(), &energyValues[0], &acceptanceValues_psd_charge_cut[0]);
    TGraph gr_acceptance_all_cut(energyValues.size(), &energyValues[0], &acceptanceValues_all_cut[0]);

    gr_acceptance_gometric_cut.SetName("gr_acceptance_gometric_cut");
    gr_acceptance_maxElateral_cut.SetName("gr_acceptance_maxElateral_cut");
    gr_acceptance_maxBarLayer_cut.SetName("gr_acceptance_maxBarLayer_cut");
    gr_acceptance_BGOTrackContainment_cut.SetName("gr_acceptance_BGOTrackContainment_cut");
    gr_acceptance_nBarLayer13_cut.SetName("gr_acceptance_nBarLayer13_cut");
    gr_acceptance_maxRms_cut.SetName("gr_acceptance_maxRms_cut");
    gr_acceptance_track_selection_cut.SetName("gr_acceptance_track_selection_cut");
    gr_acceptance_xtrl_cut.SetName("gr_acceptance_xtrl_cut");
    gr_acceptance_psd_charge_cut.SetName("gr_acceptance_psd_charge_cut");
    gr_acceptance_all_cut.SetName("gr_acceptance_all_cut");

    gr_acceptance_gometric_cut.SetTitle("Acceptance - geometric cut");
    gr_acceptance_maxElateral_cut.SetTitle("Acceptance - maxElateral cut");
    gr_acceptance_maxBarLayer_cut.SetTitle("Acceptance - maxBarLayer cut");
    gr_acceptance_BGOTrackContainment_cut.SetTitle("Acceptance - BGOTrackContainment cut");
    gr_acceptance_nBarLayer13_cut.SetTitle("Acceptance - nBarLayer13 cut");
    gr_acceptance_maxRms_cut.SetTitle("Acceptance - maxRms cut");
    gr_acceptance_track_selection_cut.SetTitle("Acceptance - track selection cut");
    gr_acceptance_xtrl_cut.SetTitle("Acceptance - XTRL cut");
    gr_acceptance_psd_charge_cut.SetTitle("Acceptance - PSD charge selection cut");
    gr_acceptance_all_cut.SetTitle("Acceptance - all cut");


    // Write histos to file
    h_incoming.Write();
    h_gometric_cut.Write();
    h_maxElayer_cut.Write();
    h_maxBarLayer_cut.Write();
    h_BGOTrackContainment_cut.Write();
    h_nBarLayer13_cut.Write();
    h_maxRms_cut.Write();
    h_track_selection_cut.Write();
    h_xtrl_cut.Write();
    h_psd_charge_cut.Write();
    h_all_cut.Write();

    // Create output acceptance dir in the output TFile
    auto acceptanceDIr = outFile.mkdir("Acceptance");
    acceptanceDIr->cd();
    
    // Write final TGraphs
    gr_acceptance_gometric_cut.Write();
    gr_acceptance_maxElateral_cut.Write();
    gr_acceptance_maxBarLayer_cut.Write();
    gr_acceptance_BGOTrackContainment_cut.Write();
    gr_acceptance_nBarLayer13_cut.Write();
    gr_acceptance_maxRms_cut.Write();
    gr_acceptance_track_selection_cut.Write();
    gr_acceptance_xtrl_cut.Write();
    gr_acceptance_psd_charge_cut.Write();
    gr_acceptance_all_cut.Write();

    // Return to main TFile directory
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