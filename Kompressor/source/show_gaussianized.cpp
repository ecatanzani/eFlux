#include "gaussianize.h"

#include <map>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>

#include "TF1.h"

void showGaussianizedTMVAvars(
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

    auto lambda_values = _lambda_config->GetLambdaStruct();
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

    auto h_rmsLayer_gauss = GetRMSLayerHistos(energy_nbins, lambda_values, DAMPE_bgo_nLayers);

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {   
        //std::cout << "\nBin idx: " << bin_idx;
        auto bin_filter = [&bin_idx](int energy_bin) -> bool { return energy_bin == bin_idx; };
        _data_fr.Filter(bin_filter, {"energy_bin"})
                .Foreach([&h_rmsLayer_gauss, &bin_idx, &lambda_values](const std::map<double, std::vector<double>> rmslayer_gauss, const double energy_w)
                {
                    if (rmslayer_gauss.size()!=lambda_values.num+1) exit(100);
                    else
                    {
                        int lambda_idx = 0;
                        for (const auto &map_elm : rmslayer_gauss)
                        {
                            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                            {
                                //std::cout << "\nBin idx: " << bin_idx << "\t lambda_idx: " << lambda_idx << "\t layer: " << ly;
                                h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->Fill(map_elm.second[ly], energy_w);
                            }
                            ++lambda_idx;
                        }
                    }
                }, {"rmsLayer_gauss", "simu_energy_w_corr"});
    }

    // Compute the goodness
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        std::vector<std::vector<double>> goodness_rms_layer (DAMPE_bgo_nLayers, std::vector<double> (lambda_values.num +1, 999));
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            for (int lambda_idx=0; lambda_idx<=lambda_values.num; ++lambda_idx)
            {
                if (h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->GetEntries()>5)
                {
                    auto xmin = h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->GetXaxis()->GetXmin();
                    auto xmax = h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->GetXaxis()->GetXmax();
                    std::unique_ptr<TF1> fitfunc = std::make_unique<TF1>("fitfunc", "gaus", xmin, xmax);
                    h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->Fit("fitfunc", "Q");
                    auto chi2 = fitfunc->GetChisquare();
                    auto dof = fitfunc->GetNDF();
                    double skew = h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->GetSkewness();
                    if ((h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->GetRMS()/skew)>0.01 && 
                        (h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->GetRMS()/skew)<100.0 && 
                        skew<1.0 && dof)
                            goodness_rms_layer[ly][lambda_idx] = abs(chi2/dof -1);    
                }
            }
            
            std::vector<double>::iterator result = std::min_element(goodness_rms_layer[ly].begin(), goodness_rms_layer[ly].end());
            auto best_goodness_idx = std::distance(goodness_rms_layer[ly].begin(), result);
            auto best_lambda = lambda_values.start + lambda_values.step*best_goodness_idx;
            auto best_goodness = goodness_rms_layer[ly][best_goodness_idx];
            auto statistics = h_rmsLayer_gauss[bin_idx-1][best_goodness_idx][ly]->GetEntries();
            std::cout << "\n[ENERGY BIN " << bin_idx << "]\t - BGO Layer " << ly+1 << "\t - best lambda value: " << best_lambda << "\t - goodness: " << best_goodness << "\t Stat (#events) [" << statistics << "]";

        }
    }

    TFile* output_file = TFile::Open(outputPath.c_str(), "RECREATE");
    if (output_file->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outputPath << "]\n\n";
        exit(100);
    }

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        output_file->mkdir((std::string("energybin_") + std::to_string(bin_idx)).c_str());
        output_file->cd((std::string("energybin_") + std::to_string(bin_idx)).c_str());
        for (int lambda_idx=0; lambda_idx<=lambda_values.num; ++lambda_idx)
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->Write();
    }

    output_file->Close();

    if (_VERBOSE)
        std::cout << "\n\nOutput file has been written... [" << outputPath << "]\n";
}

std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetAutoRMSLayerHistos(
    const int energy_nbins, 
    const lambdas lambda_values,
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
    const lambdas lambda_values,
    const int DAMPE_bgo_nLayers)
{
    int nbins = 200;
    double hmin = 0;
    double hmax = 40;
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