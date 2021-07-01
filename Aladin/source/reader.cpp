#include "reader.h"
#include "fit.h"
#include "list_parser.h"
#include "gaussianize.h"
#include "loglikelihood.h"

#include "config.h"
#include "energy_config.h"
#include "lambda_config.h"

void reader(in_args input_args)
{
    // Parse input file list
    bool gaussianized = input_args.gaussianize ? false : true;
    std::shared_ptr<parser> evt_parser = std::make_shared<parser>(input_args.input_list, input_args.mc_flag, input_args.verbose, gaussianized);
    // Parse 'Collector' config file
    std::shared_ptr<config> _config = std::make_shared<config>(input_args.wd, input_args.mc_flag);
    // Parse local config file
    std::shared_ptr<energy_config> _energy_config = std::make_shared<energy_config>(input_args.wd);
    // Parse local lambda config file
    std::shared_ptr<lambda_config> _lambda_config = std::make_shared<lambda_config>(input_args.wd);
    
    if (input_args.gaussianize)
        gaussianize(
            evt_parser->GetEvtTree(),
            _config,
            _energy_config,
            _lambda_config,
            evt_parser->GetEvtTree()->GetEntries(),
            input_args.output_path,
            input_args.regularize_tree_path, 
            input_args.verbose,
            input_args.threads,
            input_args.mc_flag);

    else if (input_args.loglikelihood)
        buildLogLikelihoodProfile(
            evt_parser->GetEvtTree(),
            _config,
            _energy_config,
            _lambda_config,
            evt_parser->GetEvtTree()->GetEntries(),
            input_args.energybin,
            input_args.output_path, 
            input_args.regularize_tree_path,
            input_args.verbose,
            input_args.threads,
            input_args.mc_flag);

    else if (input_args.fit)
        fit(evt_parser->GetEvtTree(),
            _config,
            _energy_config,
            _lambda_config,
            evt_parser->GetEvtTree()->GetEntries(),
            input_args.energybin,
            input_args.output_path, 
            input_args.regularize_tree_path,
            input_args.verbose,
            input_args.threads,
            input_args.mc_flag);

}