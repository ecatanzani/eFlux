#include "split.h"
#include "config.h"
#include "list_parser.h"
#include "energy_config.h"

#include <memory>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>

void split(in_args input_args)
{
    std::unique_ptr<parser> evt_parser = std::make_unique<parser>(input_args.input_list, input_args.mc_flag, input_args.verbose, true);
    std::unique_ptr<config> _config = std::make_unique<config>(input_args.wd, input_args.mc_flag);
    std::unique_ptr<energy_config> _energy_config = std::make_unique<energy_config>(input_args.wd);

    if (input_args.verbose)
    {
        _config->PrintActiveFilters();
        _energy_config->PrintActiveFilters();
        std::cout << "\nTotal number of events: " << evt_parser->GetEvtTree()->GetEntries();
    }

    ROOT::EnableImplicitMT(input_args.threads);
    ROOT::RDataFrame _data_fr(*evt_parser->GetEvtTree());

    if (input_args.verbose) std::cout << "\n\nAnalysis running...\n";

    auto energy_nbins = (int)(_energy_config->GetEnergyBinning().size()) - 1;

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        auto bin_filter = [bin_idx](int energy_bin) -> bool { return energy_bin == bin_idx; };
        auto fname = input_args.output_path + std::string("energybin_") + std::to_string(bin_idx) + std::string(".root");
        _data_fr.Filter(bin_filter, {"energy_bin"}).Snapshot(evt_parser->GetEvtTree()->GetName(), fname.c_str());
    }
}