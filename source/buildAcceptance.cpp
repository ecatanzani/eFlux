/*
    Consider this repo as a manual: https://github.com/DAMPEEU/DmpTools/blob/master/HE_skimmer/code/app/main.cc#L318
*/

#include "myHeader.h"
#include "acceptance.h"
#include "energyMatch.h"

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

inline double wtsydp(
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
    TFile &outFile)
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
    load_acceptance_struct(acceptance_cuts);

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

bool maxElater_cut(
    std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf &acceptance_cuts,
    const double bgoTotalE)
{
    // Get energy maximum along X and Y views
    double ELayer_max_XZ = 0;
    double ELayer_max_YZ = 0;

    for (int lIdx = 1; lIdx < DAMPE_bgo_nLayers; lIdx += 2)
    {
        auto lEgy = bgorec->GetELayer(lIdx);
        if (lEgy > ELayer_max_XZ)
            ELayer_max_XZ = lEgy;
    }

    for (int lIdx = 1; lIdx < DAMPE_bgo_nLayers; lIdx += 1)
    {
        auto lEgy = bgorec->GetELayer(lIdx);
        if (lEgy > ELayer_max_YZ)
            ELayer_max_YZ = lEgy;
    }

    bool passed_maxELayerTotalE_cut = true;
    double MaxELayer;
    if (ELayer_max_XZ > ELayer_max_YZ)
        MaxELayer = ELayer_max_XZ;
    else
        MaxELayer = ELayer_max_YZ;
    double rMaxELayerTotalE = MaxELayer / bgoTotalE;
    if (rMaxELayerTotalE > acceptance_cuts.energy_lRatio)
        passed_maxELayerTotalE_cut = false;

    return passed_maxELayerTotalE_cut;
}

bool maxBarLayer_cut(
    std::shared_ptr<DmpEvtBgoHits> bgohits,
    const int nBgoHits)
{
    bool passed_maxBarLayer_cut = true;
    std::vector<short> barNumberMaxEBarLay1_2_3(3, -1); // Bar number of maxE bar in layer 1, 2, 3
    std::vector<double> MaxEBarLay1_2_3(3, 0);          // E of maxE bar in layer 1, 2, 3

    for (int ihit = 0; ihit < nBgoHits; ++ihit)
    {
        auto hitE = (bgohits->fEnergy)[ihit];
        auto lay = bgohits->GetLayerID(ihit);
        if (lay == 1 || lay == 2 || lay == 3)
        {
            if (hitE > MaxEBarLay1_2_3[lay - 1])
            {
                auto iBar = ((bgohits->fGlobalBarID)[ihit] >> 6) & 0x1f;
                MaxEBarLay1_2_3[lay - 1] = hitE;
                barNumberMaxEBarLay1_2_3[lay - 1] = iBar;
            }
        }
    }

    for (int j = 0; j < 3; ++j)
        if (barNumberMaxEBarLay1_2_3[j] <= 0 || barNumberMaxEBarLay1_2_3[j] == 21)
            passed_maxBarLayer_cut = false;

    return passed_maxBarLayer_cut;
}

bool BGOTrackContainment_cut(
    std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf &acceptance_cuts,
    bool passEvent)
{
    bool passed_bgo_containment_cut = false;
    double BGO_TopZ = 46;
    double BGO_BottomZ = 448;
    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);
    bgoRec_slope[1] = bgorec->GetSlopeXZ();
    bgoRec_slope[0] = bgorec->GetSlopeYZ();
    bgoRec_intercept[1] = bgorec->GetInterceptXZ();
    bgoRec_intercept[0] = bgorec->GetInterceptYZ();

    if ((bgoRec_slope[1] == 0 && bgoRec_intercept[1] == 0) ||
        (bgoRec_slope[0] == 0 && bgoRec_intercept[0] == 0))
        passEvent = false;

    TVector3 bgoRecEntrance;
    TVector3 bgoRecExit;

    double topX = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];
    double topY = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
    double bottomX = bgoRec_slope[1] * BGO_BottomZ + bgoRec_intercept[1];
    double bottomY = bgoRec_slope[0] * BGO_BottomZ + bgoRec_intercept[0];

    if (fabs(topX) < acceptance_cuts.shower_axis_delta && fabs(topY) < acceptance_cuts.shower_axis_delta && fabs(bottomX) < acceptance_cuts.shower_axis_delta && fabs(bottomY) < acceptance_cuts.shower_axis_delta)
        passed_bgo_containment_cut = true;

    return passed_bgo_containment_cut;
}