#include "set.h"
#include "list_parser.h"
#include "energy_config.h"

#include <iostream>
#include <memory>

#include "TChain.h"
#include "TRandom.h"
#include <ROOT/RDataFrame.hxx>

void buildSet(in_args input_pars) {
    std::shared_ptr<parser> list_parser = std::make_shared<parser>(input_pars.input_list, input_pars.verbose);
    std::shared_ptr<energy_config> dataset_energy_range = std::make_shared<energy_config>(input_pars.working_dir);
    
    ROOT::EnableImplicitMT(input_pars.threads);
    ROOT::RDataFrame _data_fr(*(list_parser->GetEvtTree()));

    auto set_filter_min_energy = dataset_energy_range->GetMinEvtEnergy();
    auto set_filter_max_energy = dataset_energy_range->GetMaxEvtEnergy();
    auto energy_filter = [&set_filter_min_energy, &set_filter_max_energy](const double energy) -> bool {
        auto status = false;
        auto _gev = 0.001;
        if (energy*_gev >= set_filter_min_energy && energy*_gev <= set_filter_max_energy)
            status = true;
        else
            status = false;
        return status; 
    };
    auto _fr_preselected = input_pars.split ? 
        _data_fr.Filter("evtfilter_all_cut==true").Filter(energy_filter, {"energy_corr"}).Define("tt_assign", []() ->double {return gRandom->Uniform();}) :
        _data_fr.Filter("evtfilter_all_cut==true").Filter(energy_filter, {"energy_corr"});

    if (input_pars.verbose) {
        std::cout << "\n\n**** Filter statistics ****\n";
        std::cout << "***************************\n";
        std::cout << "\nPreselected events: " << *(_fr_preselected.Count());
        std::cout << "\n\n********************";
    }

    if (input_pars.split) {
        auto _fr_train = _fr_preselected.Filter("tt_assign < 0.5");
        auto _fr_test = _fr_preselected.Filter("tt_assign >= 0.5");

        _fr_train.Snapshot(
            (list_parser->GetEvtTree())->GetName(), 
            (input_pars.output_path.substr(0, input_pars.output_path.find(".root")) + std::string("_trainset.root")).c_str());

        _fr_test.Snapshot(
            (list_parser->GetEvtTree())->GetName(), 
            (input_pars.output_path.substr(0, input_pars.output_path.find(".root")) + std::string("_testset.root")).c_str());
        
        if (input_pars.verbose) {
            std::cout << "\n\nOutput file has been written... [" << input_pars.output_path.substr(0, input_pars.output_path.find(".root")) + std::string("_trainset.root") << "]";
            std::cout << "\n\nOutput file has been written... [" << input_pars.output_path.substr(0, input_pars.output_path.find(".root")) + std::string("_testset.root") << "]\n";
        }

    }
    else {
        _fr_preselected.Snapshot((list_parser->GetEvtTree())->GetName(), input_pars.output_path);

        if (input_pars.verbose) std::cout << "\n\nOutput file has been written... [" << input_pars.output_path << "]\n";
    }
}