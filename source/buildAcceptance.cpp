/*
    Consider this repo as a manual: https://github.com/DAMPEEU/DmpTools/blob/master/HE_skimmer/code/app/main.cc#L318
*/

#include "myHeader.h"
#include "acceptance.h"

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

std::shared_ptr<DmpChain> aggregateEventsDmpChain(const std::string accInputPath, const bool verbose)
{
    // ****** Access data using DAMPE Chain ******

    // Create DmpChain object
    std::shared_ptr<DmpChain> dmpch = std::make_shared<DmpChain>("CollectionTree");

    // Add MC file list to DmpChain
    dmpch->AddFromList(getListPath(accInputPath, true).c_str());
    if (verbose)
        dmpch->GetListOfFiles()->Print();

    return dmpch;
}

std::shared_ptr<TChain> aggregateEventsTChain(const std::string accInputPath, const bool verbose)
{
    // ****** Access data using ROOT TChain ******

    // Create TChain object
    //TChain* dmpch = new TChain("CollectionTree");
    std::shared_ptr<TChain> dmpch = std::make_shared<TChain>("CollectionTree");
    //std::shared_ptr < TChain >  dmpch( new TChain("CollectionTree") );

    // Reading list of MC files
    std::ifstream input_file(getListPath(accInputPath, true).c_str());
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
            std::cout << "\nAdding " << tmp_str << " to the chain ..." << std::endl;
    }

    return dmpch;
}

inline int binarySearch(const std::vector<float> &logEBins, int l, int r, const double energy)
{
    /*
        ****** BINARY SEARCH ******
    */

    if (r >= l)
    {
        int mid = l + (r - l) / 2;
        if (energy > logEBins[mid] && energy < logEBins[mid + 1])
            return mid;
        if (logEBins[mid] > energy)
            return binarySearch(logEBins, l, mid - 1, energy);
        return binarySearch(logEBins, mid + 1, r, energy);
    }
    return -1;
}

inline int linearSearch(const std::vector<float> &logEBins, const double energy)
{
    /*
        ****** ALL VECTOR SCAN ******
    */

    unsigned int lIdx = 0;
    bool found_interval = false;
    for (unsigned int vIdx = 0; vIdx < logEBins.size() - 1; ++vIdx)
        if (energy > logEBins[vIdx] && energy < logEBins[vIdx + 1])
        {
            lIdx = vIdx;
            found_interval = true;
            break;
        }
    //std::cout << "\nEnergy: " << energy << "\tLow: " << logEBins[idx] << "\tHigh: " << logEBins[idx+1];
    if (found_interval)
        return lIdx;
    else
        return -1;
}

inline void allocateParticleEnergy(
    std::vector<unsigned int> &gFactorCounts,
    const std::vector<float> &logEBins,
    const double bgoTotalE)
{
    const double energy = bgoTotalE / 1000;
    auto idx = binarySearch(logEBins, 0, logEBins.size() - 1, energy); //More efficient vector search respect to the linear one
    //auto idx = linearSearch(logEBins,energy);
    if (idx != -1)
    {
        /*
        if ( !(energy > logEBins[idx] && energy < logEBins[idx+1]) ) 
            std:: cout << "\n Idx: " << idx << "\tLow(" << logEBins[idx] << "): " << idx << "\tHigh(" << logEBins[idx+1] << "): " << idx+1 << "\tEgy: " << energy;
        */
        ++gFactorCounts[idx];
    }
    else
        std::cerr << "\nWARNING: Could not find energy bin for current event: " << energy << " GeV";

    /*
        USE EXPECTIONS FOR THAT
    */
}

inline double wtsydp(Double_t minene, Double_t maxene, Double_t index)
{
    double dene = maxene - minene;
    return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * (dene))), 1. / (index));
}

void buildAcceptance(
    const std::string accInputPath,
    const bool verbose,
    const std::vector<float> &logEBins)
{

    gSystem->Load("libDmpEvent.so");

    //auto dmpch = aggregateEventsDmpChain(accInputPath,verbose);
    auto dmpch = aggregateEventsTChain(accInputPath, verbose);

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

    // Particle acceptance counter
    unsigned int accEvtCounter = 0;

    // Initialize the particle counter for each energy bin
    std::vector<unsigned int> gFactorCounts(logEBins.size(), 0);

    // Create and load acceptance events cuts from config file
    acceptance_conf acceptance_cuts;
    load_acceptance_struct(acceptance_cuts);

    for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx)
    {
        // Get chain event
        dmpch->GetEvent(evIdx);

        // GOOD EVENT variable
        bool passEvent = true;

        // Get event total energy
        double bgoTotalE = bgorec->GetTotalEnergy(); // Energy in GeV
        // Don't process events that didn't hit the detector
        if (bgoTotalE < (acceptance_cuts.event_energy * 1000))
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
            for (unsigned int ind = 0; ind < layerBarNumber[lay].size(); ++ind)
            {
                int ihit = layerBarIndex[lay][ind];
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
                    for (unsigned int ind = 0; ind < layerBarNumber[lay].size(); ++ind)
                    {
                        auto ihit = layerBarIndex[lay][ind];
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
                for (unsigned int ind = 0; ind < layerBarNumber[lay].size(); ++ind)
                {
                    auto ihit = layerBarIndex[lay][ind];
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
        auto Xtr = pow(sumRms, 4) * fracLayer[13] / 8000000.;

        /* ********************************* */

        if (maxElater_cut(bgorec, acceptance_cuts, bgoTotalE))
            if (maxBarLayer_cut(bgohits, nBgoHits))
                if (BGOTrackContainment_cut(bgorec, acceptance_cuts, passEvent))
                {
                    ++accEvtCounter;
                    allocateParticleEnergy(gFactorCounts, logEBins, bgoTotalE);
                }
    }

    if (verbose)
        std::cout << "\n\nFiltered events: " << accEvtCounter << "/" << nevents << "\n\n";

    exit(123);
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