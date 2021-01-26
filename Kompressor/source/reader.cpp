#include "reader.h"
#include "tmvaset.h"
#include "list_parser.h"

void reader(in_args input_args)
{
    // Parse input file list
    std::unique_ptr<parser> evt_parser = std::make_unique<parser>(
        input_args.input_list, 
        input_args.mc_flag, 
        input_args._VERBOSE);
    // Parse 'Collector' config file
    std::shared_ptr<config> _config = std::make_shared<config>(
        input_args.wd, 
        input_args.mc_flag);
    // Parse local config file
    std::shared_ptr<energy_config> _energy_config = std::make_shared<energy_config>(input_args.wd);
    // Get chain entries
    const double _entries = evt_parser->GetEvtTree()->GetEntries();
    if (input_args._VERBOSE)
    {
        _config->PrintActiveFilters();
        _energy_config->PrintActiveFilters();
        std::cout << "Total number of events: " << _entries;
    }
    if (input_args.mc_flag)
    {
        if (input_args.tmva_set)
        {
            createTMVAset(
                evt_parser->GetEvtTree(),
                _config,
                _energy_config,
                input_args.fit_tree_path,
                input_args.signal,
                input_args.output_path,
                input_args._VERBOSE,
                input_args.threads);
        }
        else
        {
            mc_reader(
                evt_parser->GetEvtTree(), 
                _config,
                _energy_config,
                input_args.fit_tree_path,
                _entries, 
                input_args.output_path, 
                input_args._VERBOSE,
                input_args.threads);
        }
    }
}