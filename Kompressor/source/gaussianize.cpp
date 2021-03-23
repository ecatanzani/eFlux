#include "gaussianize.h"

#include <vector>
#include <cmath>
#include <map>

#include "DAMPE_geo_structure.h"

#include <ROOT/RDataFrame.hxx>

void gaussianizeTMVAvars(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    std::shared_ptr<lambda_config> _lambda_config,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE,
    const unsigned int threads,
    const bool _mc)
{
    // Enable multithreading
    ROOT::EnableImplicitMT(threads);
    // Create RDF
    ROOT::RDataFrame _data_fr(*evtch);
    // Extract the energy binning
    auto energy_binning = _energy_config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;
    // Create the RDFs
    double _gev = 0.001;
    auto GetEnergyBin = [=](double energy) -> int 
    { 
        int bin_idx=0;
        for (; bin_idx<energy_nbins-1; ++bin_idx)
            if (energy * _gev >= energy_binning[bin_idx] && energy * _gev < energy_binning[bin_idx+1])
                break;
        return bin_idx+1; 
    };
    // Create the RDF with energy bin and energy weight
    auto _fr_bin_patch = _data_fr.Define("energy_bin", GetEnergyBin, {"energy_corr"});
    if (_mc)
        _fr_bin_patch = _fr_bin_patch.Define("simu_energy_w_corr", [&energy_binning](const double simu_energy_w) -> double { return simu_energy_w * pow(energy_binning[0], 2); }, {"simu_energy_w"});
    else
        _fr_bin_patch = _fr_bin_patch.Define("simu_energy_w_corr", "1.");
    // Filtering events
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

    // Gaussianize TMVA variables
    auto lambda_values = _lambda_config->GetLambdaStruct();
    if (_VERBOSE)
        _lambda_config->PrintLambdaSettings();
    
    if (_VERBOSE)
        std::cout << "\nAnalysis running...\n";

    auto gaussianize_rmslayer = [&lambda_values](const std::vector<double> input_rmslayer) -> std::map<double, std::vector<double>>
    {
        auto gaussianize_elm = [](const std::vector<double> elm, const double lambda) -> std::vector<double> 
        {
            std::vector<double> elm_cp = elm;
            for (unsigned int idx=0; idx<elm.size(); ++idx)
                elm_cp[idx] = lambda ? (exp(lambda*elm[idx])-1)/lambda : elm[idx];
            return elm_cp;
        };

        std::map<double, std::vector<double>> rmsLayer_gauss;
        for (int lambda_idx=0; lambda_idx<=lambda_values.num; ++lambda_idx)
        {
            double lambda = lambda_values.start + lambda_idx*lambda_values.step;
            rmsLayer_gauss.insert(std::pair<double, std::vector<double>>(lambda, gaussianize_elm(input_rmslayer, lambda)));
        }
        return rmsLayer_gauss;
    };

    auto _fr_preselected_gauss = _fr_preselected.Define("rmsLayer_gauss", gaussianize_rmslayer, {"rmsLayer"});

    _fr_preselected_gauss.Snapshot((std::string(evtch->GetName()) + std::string("_gauss")).c_str(), outputPath);
    
    if (_VERBOSE)
        std::cout << "\nOutput TFile has been written [" << outputPath << "]\n";
}