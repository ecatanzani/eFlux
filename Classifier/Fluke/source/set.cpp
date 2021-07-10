#include "set.h"
#include "config.h"
#include "list_parser.h"
#include "energy_config.h"

#include <memory>
#include <vector>

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <ROOT/RDataFrame.hxx>

struct lambda_values 
{
    std::vector<std::vector<double>> best_rms_lambda;
    std::vector<double> best_sumrms_lambda;
    std::vector<std::vector<double>> best_fraclayer_lambda;
    std::vector<double> best_fraclast_lambda;
    std::vector<double> best_xtrl_lambda;
};

inline void loadLambdas(
    lambda_values &lambda, 
    const int energy_nbins, 
    const std::string lambda_tree_path, 
    const bool verbose)
{
    lambda.best_rms_lambda.resize(energy_nbins);
    lambda.best_sumrms_lambda.resize(energy_nbins);
    lambda.best_fraclayer_lambda.resize(energy_nbins);
    lambda.best_fraclast_lambda.resize(energy_nbins);
    lambda.best_xtrl_lambda.resize(energy_nbins);

    TFile *lambda_file = TFile::Open(lambda_tree_path.c_str(), "READ");
    if (!lambda_file->IsOpen())
    {
        std::cerr << "\n\nError reading lambda TTree [" << lambda_tree_path << "]\n\n";
        exit(100);
    }
    else
    {
        if (verbose)
            std::cout << "\nReading lambda TTree [" << lambda_tree_path << "]\n";
    }

    TTreeReader myReader("corrections_tree", lambda_file);
    TTreeReaderValue<std::vector<double>> rms_lambda(myReader, "best_rms_lambda");
    TTreeReaderValue<double> sumrms_lambda(myReader, "best_sumrms_lambda");
    TTreeReaderValue<std::vector<double>> fraclayer_lambda(myReader, "best_fraclayer_lambda");
    TTreeReaderValue<double> fraclast_lambda(myReader, "best_fraclast_lambda");
    TTreeReaderValue<double> xtrl_lambda(myReader, "best_xtrl_lambda");

    unsigned int bin_idx = 0;
    while (myReader.Next())
    {
        lambda.best_rms_lambda[bin_idx] = *rms_lambda;
        lambda.best_sumrms_lambda[bin_idx] = *sumrms_lambda;
        lambda.best_fraclayer_lambda[bin_idx] = *fraclayer_lambda;
        lambda.best_fraclast_lambda[bin_idx] = *fraclast_lambda;
        lambda.best_xtrl_lambda[bin_idx] = *xtrl_lambda;
        ++bin_idx;
    }

    lambda_file->Close();
}

void createSet(const in_args input_args)
{
    std::shared_ptr<parser> evt_parser = std::make_unique<parser>(input_args.input_list);
    std::shared_ptr<config> _config = std::make_shared<config>(input_args.wd, input_args.mc_flag);
    std::shared_ptr<energy_config> _energy_config = std::make_shared<energy_config>(input_args.wd);

    auto energy_binning = _energy_config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;

    // Read lambda parameters
    auto lambdas = lambda_values(); 
    loadLambdas(lambdas, energy_nbins, input_args.reg_fit_path, input_args.verbose);

    // Build RDF
    ROOT::EnableImplicitMT(threads);
    ROOT::RDataFrame _data_fr(*(evt_parser->GetEvtTree()));


    double _gev = 0.001;
    auto GetEnergyBin = [=](double energy) -> int 
    { 
        int bin_idx=0;
        for (; bin_idx<energy_nbins-1; ++bin_idx)
            if (energy * _gev >= energy_binning[bin_idx] && energy * _gev < energy_binning[bin_idx+1])
                break;
        return bin_idx+1; 
    };

    auto _fr_bin_patch = _data_fr.Define("energy_bin", GetEnergyBin, {"energy_corr"});
    _fr_bin_patch = _mc ? _fr_bin_patch.Define("simu_energy_w_corr", [&energy_binning](const double simu_energy_w) -> double { return simu_energy_w * pow(energy_binning[0], 2); }, {"simu_energy_w"}) : _fr_bin_patch.Define("simu_energy_w_corr", "1.");

    auto set_filter_min_energy = _energy_config->GetSetMinEvtEnergy();
    auto set_filter_max_energy = _energy_config->GetSetMaxEvtEnergy();
    auto _fr_preselected = _fr_bin_patch.Filter("evtfilter_all_cut==true")
                                        .Filter([&set_filter_min_energy, &set_filter_max_energy, &_gev](const double energy) -> bool {
                                                        auto status = false;
                                                        if (energy*_gev >= set_filter_min_energy && energy*_gev <= set_filter_max_energy)
                                                            status = true;
                                                        else
                                                            status = false;
                                                        return status; }, {"energy_corr"});

    std::cout << "\n\n**** Filter statistics ****\n";
    std::cout << "***************************\n";
    std::cout << "\nPreselected events: " << *(_fr_preselected.Count());
    std::cout << "\n\n********************";

}