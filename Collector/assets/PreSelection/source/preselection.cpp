#include "bgo.h"
#include "chain.h"
#include "histos.h"
#include "preselection.h"
#include "energy_config.h"

#include <iostream>
#include <memory>

#include "DmpEvtHeader.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"

#include "TChain.h"

void preselection(const in_pars &input_pars) {
    
    auto evtch = GetFileChain(input_pars.input_path, input_pars.verbose);
    std::shared_ptr<energy_config> econfig = std::make_shared<energy_config>(input_pars.config_wd);

    // Register Header container
    std::shared_ptr<DmpEvtHeader> evt_header = std::make_shared<DmpEvtHeader>();
    evtch->SetBranchAddress("EventHeader", &evt_header);

    // Register BGO container
    std::shared_ptr<DmpEvtBgoHits> bgohits = std::make_shared<DmpEvtBgoHits>();
    evtch->SetBranchAddress("DmpEvtBgoHits", &bgohits);

    // Register BGO REC container
    std::shared_ptr<DmpEvtBgoRec> bgorec = std::make_shared<DmpEvtBgoRec>();
    evtch->SetBranchAddress("DmpEvtBgoRec", &bgorec);
    
    auto min_evt_energy = econfig->GetMinEvtEnergy();
    auto max_evt_energy = econfig->GetMaxEvtEnergy();

    double evt_corr_energy_gev {0};
    double evt_energy_gev {0};
    double gev {0.001};
    unsigned int perc {10};

    auto nevents {evtch->GetEntries()};

    std::shared_ptr<histos> ps_histos = std::make_shared<histos>(econfig);

    if (input_pars.verbose) std::cout << "\n\nNumber of events: " << nevents << std::endl;

    for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx) {

        evtch->GetEvent(evIdx);
        if (input_pars.verbose && !((evIdx+1)%(nevents*perc/100))) {
            std::cout << "\nProcessing ... [ " << perc << "% ]";
            perc += 10;
        }
    
        evt_energy_gev = bgorec->GetTotalEnergy()*gev;
        evt_corr_energy_gev = bgorec->GetElectronEcor()*gev;

        if (evt_corr_energy_gev>=min_evt_energy && evt_corr_energy_gev<=max_evt_energy) {
            bgofiducial_distributions(bgohits, bgorec, evt_header, evt_corr_energy_gev, ps_histos);
        }

    }

    ps_histos->Write(input_pars.output_wd, input_pars.verbose);
}