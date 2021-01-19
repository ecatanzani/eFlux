#include "reader.h"
#include "tmvaset.h"

void reader(in_args input_args)
{
    std::unique_ptr<parser> evt_parser = std::make_unique<parser>(
        input_args.input_list, 
        input_args.mc_flag, 
        input_args._VERBOSE);
    std::shared_ptr<config> _config = std::make_shared<config>(
        input_args.wd, 
        input_args.mc_flag);
    const double _entries = evt_parser->GetEvtTree()->GetEntries();
    if (input_args._VERBOSE)
    {
        _config->PrintActiveFilters();
        std::cout << "Total number of events: " << _entries;
    }
    if (input_args.mc_flag)
    {
        if (input_args.tmva_set)
        {
            createTMVAset(
                evt_parser->GetEvtTree(),
                _config,
                input_args.fit_tree_path,
                input_args.signal,
                input_args.output_path,
                input_args._VERBOSE);
        }
        else
        {
            mc_reader(
                evt_parser->GetEvtTree(), 
                _config,
                input_args.fit_tree_path,
                _entries, 
                input_args.output_path, 
                input_args._VERBOSE);
        }
    }
}