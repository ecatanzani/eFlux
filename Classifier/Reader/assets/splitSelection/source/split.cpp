#include "split.h"
#include "list_parser.h"
#include "energy_config.h"

#include <memory>
#include <ROOT/RDataFrame.hxx>

void Split(const in_args &input_args) {

    std::shared_ptr<parser> evt_parser = std::make_unique<parser>(input_args.input_list, input_args.verbose);
    std::shared_ptr<energy_config> config = std::make_shared<energy_config>(input_args.energy_config_file);

    if (input_args.verbose) {
        std::cout << "\n\nTotal number of entries in input file: " << evt_parser->GetEvtTree()->GetEntries() << std::endl;
        std::cout << "\nAnalysis running..." << std::endl;
    }

    ROOT::EnableImplicitMT(input_args.threads);
    ROOT::RDataFrame _data_fr(*(evt_parser->GetEvtTree()));

    for (int bin_idx = 1; bin_idx <= static_cast<int>(config->GetEnergyBinning().size()) - 1; ++bin_idx) {
        auto bin_filter = [bin_idx](int energy_bin) -> bool { return energy_bin == bin_idx; };

        std::string fname = input_args.output_directory + std::string("/") + evt_parser->GetEvtTree()->GetName() + std::string("_energybin_") + std::to_string(bin_idx) + std::string(".root");
        auto tmp_entries = *(_data_fr.Filter(bin_filter, {"energy_bin"}).Count());
        if (tmp_entries) {
            if (input_args.verbose) std::cout << "\nOutput TFile has been written [" << fname << "]\t entries: " << tmp_entries;
            _data_fr.Filter(bin_filter, {"energy_bin"}).Snapshot(evt_parser->GetEvtTree()->GetName(), fname.c_str());
        }
    }    
    

}