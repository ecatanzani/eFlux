#include "parser.h"
#include "config.h"
#include "efficiency.h"

#include "TFile.h"
#include "TChain.h"
#include <ROOT/RDataFrame.hxx>

void BuildEfficiencyHistos(in_pars input_pars)
{
    std::unique_ptr<parser> evt_parser = std::make_unique<parser>(
        input_pars.input_path,
        input_pars.mc_flag,
        input_pars.verbose);
    std::unique_ptr<config> _config = std::make_unique<config>(
        input_pars.wd, 
        input_pars.mc_flag);
    const double _entries = evt_parser->GetEvtTree()->GetEntries();
    if (input_pars.verbose)
    {
        _config->PrintActiveFilters();
        std::cout << "Total number of events: " << _entries;
    }

    // Enable multithreading
    ROOT::EnableImplicitMT(input_pars.threads);
    // Create RDF
    ROOT::RDataFrame _data_fr(*evt_parser->GetEvtTree());
    // Extract the energy binning
    auto energy_binning = _config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;


}