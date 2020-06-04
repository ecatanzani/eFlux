#include "data_loop.h"
#include "aggregate_events.h"
#include "wtsydp.h"

TH1D evLoop(
    const std::vector<float> &logEBins,
    const std::string inputPath,
    TFile &outFile,
    const bool verbose)
{
    //auto dmpch = aggregateEventsDmpChain(inputPath,verbose);
    auto dmpch = aggregateEventsTChain(inputPath, verbose);

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

    TH1D h_trigger("h_trigger", "Energy Distribution of the triggered particles", logEBins.size() - 1, &(logEBins[0]));

    h_trigger.Sumw2();



}