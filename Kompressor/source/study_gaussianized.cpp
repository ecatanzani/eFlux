#include "gaussianize.h"

#include <map>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>

#include "TF1.h"

void studyGaussianizedTMVAvars(
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
    auto rms_lambda_values = _lambda_config->GetRMSLambdaStruct();
    auto elf_lambda_values = _lambda_config->GetELFLambdaStruct();
    if (_VERBOSE)
    {
        _lambda_config->PrintLambdaSettings();
        std::cout << "\nAnalysis running..." << std::endl;
    }

#if 0
    if (_VERBOSE)
        std::cout << "\nRading histos boundaries from data frame...";
    
    std::vector<std::vector<std::vector<double>>> rms_boundaries (energy_nbins);

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        auto bin_filter = [&bin_idx](int energy_bin) -> bool { return energy_bin == bin_idx; };
        rms_boundaries[bin_idx-1] = std::vector<std::vector<double>> (DAMPE_bgo_nLayers, std::vector<double> (2, 0));
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            rms_boundaries[bin_idx-1][ly][0] = *(_data_fr.Filter(bin_filter, {"energy_bin"})
                            .Define("rms_single_layer", [&ly](const std::vector<double> rms) -> double {return rms[ly];}, {"rmsLayer"})
                            .Min<double>("rms_single_layer"));
            rms_boundaries[bin_idx-1][ly][1] = *(_data_fr.Filter(bin_filter, {"energy_bin"})
                            .Define("rms_single_layer", [&ly](const std::vector<double> rms) -> double {return rms[ly];}, {"rmsLayer"})
                            .Max<double>("rms_single_layer"));
        }
    }
#endif

    //auto h_rmsLayer_gauss = GetRMSLayerHistos(energy_nbins, lambda_values, DAMPE_bgo_nLayers, rms_boundaries);
    auto h_rmsLayer_gauss = GetRMSLayerHistos(energy_nbins, rms_lambda_values, DAMPE_bgo_nLayers);
    auto h_fracLayer_gauss = GetFracLayerHistos(energy_nbins, elf_lambda_values, DAMPE_bgo_nLayers);

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {   
        auto bin_filter = [&bin_idx](int energy_bin) -> bool { return energy_bin == bin_idx; };
        _data_fr.Filter(bin_filter, {"energy_bin"})
                .Foreach([&h_rmsLayer_gauss, &h_fracLayer_gauss, &bin_idx] (
                    const std::map<double, std::vector<double>> rmslayer_gauss,
                    const std::map<double, std::vector<double>> fraclayer_gauss, 
                    const double energy_w)
                {
                    int lambda_idx = 0;
                    for(auto map_elm = rmslayer_gauss.begin(); map_elm!=rmslayer_gauss.end(); ++map_elm)
                    {
                        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                            h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->Fill((*map_elm).second[ly], energy_w);
                        ++lambda_idx;
                    }
                    lambda_idx = 0;
                    for(auto map_elm = fraclayer_gauss.begin(); map_elm!=fraclayer_gauss.end(); ++map_elm)
                    {
                        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                            h_fracLayer_gauss[bin_idx-1][lambda_idx][ly]->Fill((*map_elm).second[ly], energy_w);
                        ++lambda_idx;
                    }
                }, {"rmsLayer_gauss", "fracLayer_gauss", "simu_energy_w_corr"});
    }
    
    // Compute the goodness
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        std::vector<double> goodness_rms_layer (DAMPE_bgo_nLayers, 999);
        std::vector<double> goodness_frac_layer (DAMPE_bgo_nLayers, 999);

        std::vector<double> rms_best_lambda (DAMPE_bgo_nLayers, rms_lambda_values.start);
        std::vector<double> fraclayer_best_lambda (DAMPE_bgo_nLayers, elf_lambda_values.start);
        std::vector<int> rms_best_lambda_idx (DAMPE_bgo_nLayers, 0);
        std::vector<int> fraclayer_best_lambda_idx (DAMPE_bgo_nLayers, 0);
        std::vector<double> rms_best_lambda_gmean (DAMPE_bgo_nLayers, 0);
        std::vector<double> rms_best_lambda_gsigma (DAMPE_bgo_nLayers, 0);
        std::vector<double> fraclayer_best_lambda_gmean (DAMPE_bgo_nLayers, 0);
        std::vector<double> fraclayer_best_lambda_gsigma (DAMPE_bgo_nLayers, 0);

        int rms_statistics = h_rmsLayer_gauss[bin_idx-1][0][0]->GetEntries();
        int fraclast_statistics = h_fracLayer_gauss[bin_idx-1][0][0]->GetEntries();

        // Get RMS info
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            for (int lambda_idx=0; lambda_idx<=rms_lambda_values.num; ++lambda_idx)
                extract_lamda_info(
                    h_rmsLayer_gauss, 
                    goodness_rms_layer, 
                    rms_best_lambda, 
                    rms_best_lambda_idx,
                    rms_best_lambda_gmean,
                    rms_best_lambda_gsigma, 
                    rms_lambda_values.start,
                    rms_lambda_values.step, 
                    bin_idx, 
                    lambda_idx, 
                    ly);
           
        // Get fraclayer info
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            for (int lambda_idx=0; lambda_idx<=elf_lambda_values.num; ++lambda_idx)
                extract_lamda_info(
                    h_fracLayer_gauss, 
                    goodness_frac_layer, 
                    fraclayer_best_lambda, 
                    fraclayer_best_lambda_idx,
                    fraclayer_best_lambda_gmean,
                    fraclayer_best_lambda_gsigma,
                    elf_lambda_values.start,
                    elf_lambda_values.step,
                    bin_idx, 
                    lambda_idx, 
                    ly);
            
        if (_VERBOSE)
        {
            std::cout << "\n\nEnergy bin [" << bin_idx << "]\n";
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                std::cout << "\n(RMS) - BGO Layer " << ly+1 << "\t - best lambda value: " << rms_best_lambda[ly] << "\t - goodness: " << 
            goodness_rms_layer[ly] << "\t Stat (#events) [" << rms_statistics << "]" << "\t best lambda idx: " << rms_best_lambda_idx[ly];
            std::cout << std::endl;
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                std::cout << "\n(EnergyLastLayerFraction) - BGO Layer " << ly+1 << "\t - best lambda value: " << fraclayer_best_lambda[ly] << "\t - goodness: " << 
            goodness_frac_layer[ly] << "\t Stat (#events) [" << fraclast_statistics << "]" << "\t best lambda idx: " << fraclayer_best_lambda_idx[ly];
        }
    }    

    // Write output file
    TFile* output_file = TFile::Open(outputPath.c_str(), "RECREATE");
    if (output_file->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outputPath << "]\n\n";
        exit(100);
    }

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        output_file->mkdir((std::string("RMS/energybin_") + std::to_string(bin_idx)).c_str());
        output_file->cd((std::string("RMS/energybin_") + std::to_string(bin_idx)).c_str());
        for (int lambda_idx=0; lambda_idx<=rms_lambda_values.num; ++lambda_idx)
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->Write();
    
        output_file->mkdir((std::string("ELF/energybin_") + std::to_string(bin_idx)).c_str());
        output_file->cd((std::string("ELF/energybin_") + std::to_string(bin_idx)).c_str());
        for (int lambda_idx=0; lambda_idx<=elf_lambda_values.num; ++lambda_idx)
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                h_fracLayer_gauss[bin_idx-1][lambda_idx][ly]->Write();
    }

    output_file->Close();

    if (_VERBOSE)
        std::cout << "\n\nOutput file has been written... [" << outputPath << "]\n";
}

std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetAutoRMSLayerHistos(
    const int energy_nbins, 
    const rms_lambdas lambda_values,
    const int DAMPE_bgo_nLayers,
    const std::vector<std::vector<std::vector<double>>> &rms_boundaries)
{
    int nbins = 100;
    double kScale = 2;
    std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> h_rms_layer (energy_nbins);
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        h_rms_layer[bin_idx-1] = std::vector<std::vector<std::shared_ptr<TH1D>>> (lambda_values.num+1);
        auto lambda = lambda_values.start;
        auto gaussianize_elm = [&lambda](const double elm) -> double { return lambda ? (exp(lambda*elm)-1)/lambda : elm; };
        for (int lambda_idx=0; lambda_idx<=lambda_values.num; ++lambda_idx)
        {
            if (lambda_idx)
                lambda += lambda_values.step;
            h_rms_layer[bin_idx-1][lambda_idx] = std::vector<std::shared_ptr<TH1D>> (DAMPE_bgo_nLayers);
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                std::string str_lambda = lambda<0 ? std::string("neg_") + std::to_string(abs(lambda)) : std::to_string(lambda);
                std::string h_name = std::string("h_rms_energybin_") + std::to_string(bin_idx) + std::string("_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
                h_rms_layer[bin_idx-1][lambda_idx][ly] = std::make_shared<TH1D>(h_name.c_str(), h_name.c_str(), nbins, gaussianize_elm(rms_boundaries[bin_idx-1][ly][0])/kScale, gaussianize_elm(rms_boundaries[bin_idx-1][ly][1])*kScale);
            }
        }
    }
    return h_rms_layer;
}

std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetRMSLayerHistos(
    const int energy_nbins, 
    const rms_lambdas lambda_values,
    const int DAMPE_bgo_nLayers)
{
    int nbins = 200;
    double hmin = 0;
    double hmax = 100;
    std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> h_rms_layer (energy_nbins);
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        h_rms_layer[bin_idx-1] = std::vector<std::vector<std::shared_ptr<TH1D>>> (lambda_values.num+1);
        auto lambda = lambda_values.start;
        for (int lambda_idx=0; lambda_idx<=lambda_values.num; ++lambda_idx)
        {
            if (lambda_idx)
                lambda += lambda_values.step;
            h_rms_layer[bin_idx-1][lambda_idx] = std::vector<std::shared_ptr<TH1D>> (DAMPE_bgo_nLayers);
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                std::string str_lambda = lambda<0 ? std::string("neg_") + std::to_string(abs(lambda)) : std::to_string(lambda);
                std::string h_name = std::string("h_rms_energybin_") + std::to_string(bin_idx) + std::string("_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
                h_rms_layer[bin_idx-1][lambda_idx][ly] = std::make_shared<TH1D>(h_name.c_str(), h_name.c_str(), nbins, hmin, hmax);
            }
        }
    }
    return h_rms_layer;
}

std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetFracLayerHistos(
    const int energy_nbins, 
    const energylastfraction_lambdas lambda_values,
    const int DAMPE_bgo_nLayers)
{
    int nbins = 1000;
    double hmin = 0;
    double hmax = 50;
    std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> h_frac_layer (energy_nbins);
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        h_frac_layer[bin_idx-1] = std::vector<std::vector<std::shared_ptr<TH1D>>> (lambda_values.num+1);
        auto lambda = lambda_values.start;
        for (int lambda_idx=0; lambda_idx<=lambda_values.num; ++lambda_idx)
        {
            if (lambda_idx)
                lambda += lambda_values.step;
            h_frac_layer[bin_idx-1][lambda_idx] = std::vector<std::shared_ptr<TH1D>> (DAMPE_bgo_nLayers);
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                std::string str_lambda = lambda<0 ? std::string("neg_") + std::to_string(abs(lambda)) : std::to_string(lambda);
                std::string h_name = std::string("h_fraclayer_energybin_") + std::to_string(bin_idx) + std::string("_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
                h_frac_layer[bin_idx-1][lambda_idx][ly] = std::make_shared<TH1D>(h_name.c_str(), h_name.c_str(), nbins, hmin, hmax);
            }
        }
    }
    return h_frac_layer;
}