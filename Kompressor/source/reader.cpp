#include "reader.h"
#include "tmvaset.h"
#include "kompress.h"
#include "list_parser.h"
#include "gaussianize.h"

#include "config.h"
#include "energy_config.h"
#include "lambda_config.h"

void reader(in_args input_args)
{
    // Parse input file list
    std::unique_ptr<parser> evt_parser = std::make_unique<parser>(
        input_args.input_list, 
        input_args.mc_flag, 
        input_args._VERBOSE,
        input_args.show_gaussianized);
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

    if (input_args.tmva_set)
        createTMVAset(
            evt_parser->GetEvtTree(),
            _config,
            _energy_config,
            input_args.fit_tree_path,
            input_args.signal,
            input_args.output_path,
            input_args._VERBOSE,
            input_args.no_split,
            input_args.no_split_test,
            input_args.threads,
            input_args.mc_flag);
    else if (input_args.gaussianize)
    {
        std::shared_ptr<lambda_config> _lambda_config = std::make_shared<lambda_config>(input_args.wd);
        gaussianizeTMVAvars(
            evt_parser->GetEvtTree(),
            _config,
            _energy_config,
            _lambda_config,
            _entries,
            input_args.output_path, 
            input_args._VERBOSE,
            input_args.threads,
            input_args.mc_flag);
    }
    else if (input_args.show_gaussianized)
    {
        std::shared_ptr<lambda_config> _lambda_config = std::make_shared<lambda_config>(input_args.wd);
        showGaussianizedTMVAvars(
            evt_parser->GetEvtTree(),
            _config,
            _energy_config,
            _lambda_config,
            _entries,
            input_args.output_path, 
            input_args._VERBOSE,
            input_args.threads,
            input_args.mc_flag);
    }
    else
        kompress(
            evt_parser->GetEvtTree(), 
            _config,
            _energy_config,
            input_args.fit_tree_path,
            _entries, 
            input_args.output_path, 
            input_args._VERBOSE,
            input_args.threads,
            input_args.mc_flag);

}