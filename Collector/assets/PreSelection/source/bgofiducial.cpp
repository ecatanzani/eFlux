#include "bgofiducial.h"

#include "Dmp/DmpBgoContainer.h"

#include "DmpEvtHeader.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"

#include "TH1D.h"

inline bool check_trigger(const std::shared_ptr<DmpEvtHeader> evt_header) {
	auto trigger_mip1 = evt_header->GeneratedTrigger(1);
	auto trigger_mip2 = evt_header->GeneratedTrigger(2);
	auto trigger_HET = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
	auto trigger_LET = evt_header->GeneratedTrigger(4) && evt_header->EnabledTrigger(4);
	auto trigger_MIP = trigger_mip1 || trigger_mip2;
	auto trigger_general = trigger_MIP || trigger_HET || trigger_LET;
	return trigger_HET;
}

void bgofiducial_distributions(const std::string output_path, std::shared_ptr<TChain> evtch, const bool verbose) {

    // Register Header container
	std::shared_ptr<DmpEvtHeader> evt_header = std::make_shared<DmpEvtHeader>();
	evtch->SetBranchAddress("EventHeader", &evt_header);

    // Register BGO container
	std::shared_ptr<DmpEvtBgoHits> bgohits = std::make_shared<DmpEvtBgoHits>();
	evtch->SetBranchAddress("DmpEvtBgoHits", &bgohits);

    // Register BGO REC container
	std::shared_ptr<DmpEvtBgoRec> bgorec = std::make_shared<DmpEvtBgoRec>();
	evtch->SetBranchAddress("DmpEvtBgoRec", &bgorec);


    double layer_min_energy {0};           //Minimum energy per layer
    auto nevents {evtch->GetEntries()};

    std::unique_ptr<TH1D> h_energy_fraction = std::make_unique<TH1D>("h_energy_fraction", "Energy fraction; Energy fraction; counts", 100, 0, 1);
    
    unsigned int step {10000};
    if (verbose) std::cout << "\n\nNumber of events: " << nevents << std::endl;
    for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx) {
        evtch->GetEvent(evIdx);
        if (verbose && !((evIdx+1)%step))
            std::cout << "\nNumber of processed events [" << evIdx+1 << "]";
        
        if (check_trigger(evt_header)) {
            std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
            bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), layer_min_energy);
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                h_energy_fraction->Fill(elm_energy_fraction);
        }
    }

    TFile* outfile = TFile::Open(output_path.c_str(), "RECREATE");
    if (!outfile->IsOpen()) {
        std::cerr << "\n\nError writing output file [" << output_path << "]\n\n";
        exit(100);
    }
    
    h_energy_fraction->Write();

    outfile->Close();

    if (verbose) std::cout << "\n\nOutput file has been written [" << output_path << "]\n\n";
}