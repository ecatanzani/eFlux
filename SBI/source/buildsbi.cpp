#include "sbi.h"
#include "timing.h"
#include "buildsbi.h"
#include "aggregate_events.h"

#include "DmpChain.h"

void buildSBI(const in_pars input_pars)
{
    // Create output file
    TFile* outfile = TFile::Open(input_pars.output_path.c_str(), "RECREATE");
    if (!outfile->IsOpen())
    {
        std::cerr << "\n\nERROR: output file not created [" << input_pars.output_path << "]\n\n";
        exit(100);
    }

    // Aggregate files
    std::shared_ptr<DmpChain> dmpch;
	event_collector collector(input_pars.input_path, input_pars.verbose);
	if (collector.GetChainStatus())
		dmpch = collector.GetChain();
	else
	{
		std::cerr << "\n\nERROR: Corrupted DmpChain object...";
		exit(100);
	}

    // Create sbi container
    std::unique_ptr<sbi> evt_sbi = std::make_unique<sbi>();
    // Create timing container
    std::unique_ptr<timing> evt_time = std::make_unique<timing>();

    if (input_pars.verbose)
    {
        std::cout << "\n\nNumber of events: " << dmpch->GetEntries() << "\n\n";
        std::cout << "\nAnalysis running...\n\n";
    }

    for (unsigned int ev_idx=1; ev_idx<dmpch->GetEntries(); ++ev_idx)
    {
        auto pev = dmpch->GetDmpEvent(ev_idx);
        auto header = pev->pEvtHeader();
        auto evtatt = pev->pEvtAttitude();
        evt_time->SetCurrSec(header->GetSecond());
        evt_time->CheckRepeated();

        // Event check
        if (!evt_sbi->SetSBIStatus(header, evtatt, ev_idx))
            continue;

        

    }

    outfile->Close();
    outfile->Write();
}