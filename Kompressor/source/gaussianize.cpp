#include "gaussianize.h"

#include <vector>
#include <cmath>

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
    
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        auto lambda = lambda_values.start;
        for (; lambda<=lambda_values.end; lambda+=lambda_values.step)
        {
            auto gaussianize = [&lambda](const std::vector<double> x) -> std::vector<double> 
            { 
                std::vector<double> x2 = x;
                for (unsigned int idx=0; idx<x2.size(); ++idx)
                    x2[idx] = lambda ? (exp(lambda*x[idx])-1)/lambda : x[idx];
                return x2;
            };
            std::string str_lambda = lambda<0 ? std::string("neg_") + std::to_string(abs(lambda)) : std::to_string(lambda);
            std::string leaf_name = std::string("rmsLayer_energybin_") + std::to_string(bin_idx) + std::string("_lambda_") + str_lambda;
            leaf_name.erase(leaf_name.find_last_not_of('0'), std::string::npos);
            _fr_preselected = _fr_preselected.Define(leaf_name, gaussianize, {"rmsLayer"});
        }
    }    
    
    // Create histos
    if (_VERBOSE)
        std::cout << "\nGaussianizing TMVA variables...\n";

    std::vector<std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>>> h_rmsLayer_lambda (energy_nbins);

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        if (_VERBOSE)
            std::cout << "\n[INFO] Building histos - energy bin [" << bin_idx << "]";
        auto bin_filter = [&bin_idx](int energy_bin) -> bool { return energy_bin == bin_idx; };
        h_rmsLayer_lambda[bin_idx-1] = std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> (lambda_values.num);
        auto lambda = lambda_values.start;
        for (int lambda_idx=0; lambda_idx<lambda_values.num; ++lambda_idx)
        {
            if (lambda_idx)
                lambda += lambda_values.step;
            h_rmsLayer_lambda[bin_idx-1][lambda_idx] = std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers);
            std::string str_lambda = lambda<0 ? std::string("neg_") + std::to_string(abs(lambda)) : std::to_string(lambda);
            std::string leaf_name = std::string("rmsLayer_energybin_") + std::to_string(bin_idx) + std::string("_lambda_") + str_lambda;
            leaf_name.erase(leaf_name.find_last_not_of('0'), std::string::npos);
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                auto get_layer_info = [&ly](std::vector<double> rms_layer) -> double { return rms_layer[ly]; };
                h_rmsLayer_lambda[bin_idx-1][lambda_idx][ly] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                                .Define("rms_layer", get_layer_info, {leaf_name})
                                                                                .Histo1D<double, double>("rms_layer", "simu_energy_w_corr");
                h_rmsLayer_lambda[bin_idx-1][lambda_idx][ly]->SetName(
                    (std::string("h_rmsLayer_energybin_") + std::to_string(bin_idx) + std::string("_lambda_") + std::to_string(lambda) + "_layer_" + std::to_string(ly)).c_str());
                h_rmsLayer_lambda[bin_idx-1][lambda_idx][ly]->SetTitle((std::string("rmsLayer - energybin ") + std::to_string(bin_idx) + std::string(" - lambda ") + std::to_string(lambda) + " - layer " + std::to_string(ly)).c_str());
            }
        }
    }

    //ComputeGoodness(h_rms);

    if (_VERBOSE)
        std::cout << "\n\nWriting output file...\n";
    TFile* output_file = TFile::Open(outputPath.c_str(), "RECREATE");
    if (output_file->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outputPath << "]\n\n";
        exit(100);
    }

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        auto lambda_idx = 0;
        auto lambda = lambda_values.start;
        for (; lambda<=lambda_values.end; lambda+=lambda_values.step)
        {
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                h_rmsLayer_lambda[bin_idx-1][lambda_idx][ly]->Write();
            }
            ++lambda_idx;
        }
    }

    output_file->Close();

}

void ComputeGoodness(std::vector<std::shared_ptr<TH1D>> h_rms)
{
    std::vector<double> goodness_rms (h_rms.size(), 999);

}