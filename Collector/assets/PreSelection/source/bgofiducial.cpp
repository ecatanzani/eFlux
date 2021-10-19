#include "bgofiducial.h"

#include "Dmp/DmpBgoContainer.h"

#include "DmpEvtHeader.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"

#include "TH1D.h"

void bgofiducial_distributions(TFile* outfile, std::shared_ptr<TChain> evtch, const bool verbose) {

    // Register Header container
	std::shared_ptr<DmpEvtHeader> evt_header = std::make_shared<DmpEvtHeader>();
	evtch->SetBranchAddress("EventHeader", &evt_header);

    // Register BGO container
	std::shared_ptr<DmpEvtBgoHits> bgohits = std::make_shared<DmpEvtBgoHits>();
	evtch->SetBranchAddress("DmpEvtBgoHits", &bgohits);

    // Register BGO REC container
	std::shared_ptr<DmpEvtBgoRec> bgorec = std::make_shared<DmpEvtBgoRec>();
	evtch->SetBranchAddress("DmpEvtBgoRec", &bgorec);

    std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();

    double layer_min_energy {0};           //Minimum energy per layer
    auto nevents  {evtch->GetEntries()};

    std::unique_ptr<TH1D> h_energy_fraction = std::make_unique<TH1D>("h_energy_fraction", "Energy fraction; Energy fraction; counts", 100, 0, 1);

    for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx) {
        evtch->GetEvent(evIdx);
        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), layer_min_energy);

        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
            h_energy_fraction->Fill(elm_energy_fraction);
    }

    h_energy_fraction->Write();

    outfile->Close();

}