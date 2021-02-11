#include "sbi.h"
#include "utils.h"
#include "buildsbi.h"
#include "container.h"
#include "aggregate_events.h"

#include "DmpEvtNudRaw.h"
#include "DmpEvtHeader.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtPsdHits.h"
#include "DmpEvtStkReco.hh"
#include "DmpEvtAttitude.h"

#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"

void buildSBI(const in_pars input_pars)
{
    // Aggregate files
    std::shared_ptr<TChain> dmpch;
    event_collector collector(input_pars.input_path, input_pars.verbose);
    if (collector.GetChainStatus())
		dmpch = collector.GetChain();
	else
	{
		std::cerr << "\n\nERROR: Corrupted DmpChain object...";
		exit(100);
	}
    
    // Load output file
    TFile outfile(input_pars.output_path.c_str(), "NEW", "SBI Output File");
    if (!outfile.IsOpen())
    {
        std::cerr << "\n\nERROR: output file not created [" << input_pars.output_path << "]\n\n";
        exit(100);
    }

    // Register Header container
	std::shared_ptr<DmpEvtHeader> header = std::make_shared<DmpEvtHeader>();
	dmpch->SetBranchAddress("EventHeader", &header);

    // Register BGO REC container
	std::shared_ptr<DmpEvtBgoRec> bgorec = std::make_shared<DmpEvtBgoRec>();
	dmpch->SetBranchAddress("DmpEvtBgoRec", &bgorec);
  
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

    // Register NUD container
    std::shared_ptr<DmpEvtNudRaw> nudraw = std::make_shared<DmpEvtNudRaw>();
	dmpch->SetBranchAddress("DmpEvtNudRaw", &nudraw);

    // Register attitude container
	std::shared_ptr<DmpEvtAttitude> attitude = std::make_shared<DmpEvtAttitude>();
	dmpch->SetBranchAddress("EvtAttitudeContainer", &attitude);

    // Create sbi container
    std::shared_ptr<sbi> evt_sbi = std::make_shared<sbi>();
    
    // Create timing container
    std::shared_ptr<container> sbi_second_container = std::make_shared<container>();
    
    if (input_pars.verbose)
    {
        std::cout << "\n\nNumber of events: " << dmpch->GetEntries() << std::endl;
        std::cout << "\nAnalysing RawData...\n";
    }
    
    double acctime;
    int kStep = 10;
    auto nevents = dmpch->GetEntries();

    for (unsigned int ev_idx=1; ev_idx<nevents; ++ev_idx)
    {
        if (input_pars.verbose)
            UpdateProcessStatus(ev_idx, kStep, nevents);

        // Read tree event
		dmpch->GetEvent(ev_idx);
        
        // tmp second
        auto tmp_second = header->GetSecond();

        /*
            - Set current second 
            - Checks if the second cnahges respect to the previous one
            - Checks if the seconds has been repeated
        */
        sbi_second_container->evt_time.SetCurrSec(tmp_second);
        
        // Event check
        if (!sbi_second_container->SetSBIStatus(header, attitude, ev_idx))
            continue;
        
        // Update SBI second
        sbi_second_container->SetSecond();
        // Update SBI event number
        sbi_second_container->UpdateEventNumber();

        // First event
        if (ev_idx==1)
        {
            // Set lifetime for the first event
            acctime = (header->GetMillisecond() + _msdeadtime)/(double)1000;
            // Update SBI run number
            sbi_second_container->UpdateRunNumber();
        }
        // Other events
        else if (ev_idx!=1 && sbi_second_container->evt_time.CheckNewSec())
        {
            if (sbi_second_container->CheckRepeatedSeconds())
            {
                // Compute SBI livetime
                sbi_second_container->ComputeLiveTime();
                // Fill SBI tree
                evt_sbi->Fill(sbi_second_container, attitude);
                // Reset SBI clasds
                evt_sbi->Reset();
                // Cleanup SBI container
                sbi_second_container->CleanUp(ev_idx);
            }
        }
        else if (ev_idx!=1 && !sbi_second_container->evt_time.CheckNewSec())
        {
            // Set SBI last event number
            sbi_second_container->SetLastEventNumber(ev_idx);
            // Livetime
            sbi_second_container->UpdateDeadTime(header->GetMillisecond());
            // Update lifetime
            acctime += (header->GetDeltaTime())/(double)1000;
            // Update trigger
            sbi_second_container->UpdateTrigger(header);
            // Update subdetectors occupancy
            sbi_second_container->UpdateSubDetectorsOccupancy(
                stktracks,
                bgorec,
                psdhits,
                nudraw);
            // Update subdetector status
            sbi_second_container->UpdateSubDetectorStatus(header);
            // Check SAA
            sbi_second_container->CheckSAA(attitude);
        }
        
        // Set previous second
        sbi_second_container->evt_time.SetPrevSec(tmp_second);
    }

    // Last fill before closing the file
    sbi_second_container->SetSBIStatus(false);
    // Fill SBI tree
    evt_sbi->Fill(sbi_second_container, attitude);

    if (input_pars.verbose)
    {
        std::cout << "\n\nLooping done...";
        std::cout << "\nComputed lifetime: " << acctime << std::endl;
    }

    evt_sbi->Write(outfile);
    outfile.Close();
}