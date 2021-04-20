#include "gaussianize.h"

#include <map>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>

#include "TF1.h"
#include "TFile.h"

void fitGaussianizedTMVAvars(
    const std::string input_file,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    std::shared_ptr<lambda_config> _lambda_config,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE,
    const unsigned int threads,
    const bool _mc)
{
    // Extract the energy binning
    auto energy_binning = _energy_config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;
    // Extract lambda config
    auto lambda_values = _lambda_config->GetLambdaStruct();
    if (_VERBOSE) _lambda_config->PrintLambdaSettings();
    
    // Get histo structure
    auto h_rmsLayer_gauss = GetRMSLayerHistosPointers(energy_nbins, lambda_values, DAMPE_bgo_nLayers);

    // Read input file
    TFile *infile = TFile::Open(input_file.c_str(), "READ");
    if (infile->IsZombie())
    {
        std::cerr << "\n\nError opening input file [" << input_file << "]\n\n";
        exit(100);
    }

    // Read input file structure
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        auto lambda = lambda_values.start;
        for (int lambda_idx=0; lambda_idx<=lambda_values.num; ++lambda_idx)
        {
            if (lambda_idx) lambda += lambda_values.step;
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                std::string str_lambda = lambda<0 ? std::string("neg_") + std::to_string(abs(lambda)) : std::to_string(lambda);
                std::string h_name = std::string("h_rms_energybin_") + std::to_string(bin_idx) + std::string("_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
                h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly] = std::shared_ptr<TH1D>(static_cast<TH1D*>(infile->Get((std::string("energybin_") + std::to_string(bin_idx) + std::string("/") + h_name).c_str())));
                h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->SetDirectory(0);
            }
        }
    }

    infile->Close();

    // Create output log file
    std::ofstream _log("lambdafit_boxcox.log");

    // Compute the goodness
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        std::vector<double> goodness_rms_layer (DAMPE_bgo_nLayers, 999);
        std::vector<double> best_lambda (DAMPE_bgo_nLayers, lambda_values.start);
        int statistics = h_rmsLayer_gauss[bin_idx-1][0][0]->GetEntries();

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
                    {
                        double tmp_goodness = chi2/dof;
                        if (goodness_rms_layer[ly]==999)
                        {
                            goodness_rms_layer[ly] = tmp_goodness;
                            best_lambda[ly] = lambda_values.start + lambda_values.step*lambda_idx;
                        }
                        else
                        {
                            if (abs(tmp_goodness-1) < abs(goodness_rms_layer[ly]-1))
                            {
                                goodness_rms_layer[ly] = tmp_goodness;
                                best_lambda[ly] = lambda_values.start + lambda_values.step*lambda_idx;
                            }
                        }   
                    }
                }
            }
            
            std::cout << "\n[ENERGY BIN " << bin_idx << "]\t - BGO Layer " << ly+1 << "\t - best lambda value: " << best_lambda[ly] << "\t - goodness: " << 
            goodness_rms_layer[ly] << "\t Stat (#events) [" << statistics << "]";
            _log << "\n[ENERGY BIN " << bin_idx << "]\t - BGO Layer " << ly+1 << "\t - best lambda value: " << best_lambda[ly] << "\t - goodness: " << 
            goodness_rms_layer[ly] << "\t Stat (#events) [" << statistics << "]";
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

std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetRMSLayerHistosPointers(
    const int energy_nbins, 
    const lambdas lambda_values,
    const int DAMPE_bgo_nLayers)
{
    std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> h_rms_layer (energy_nbins);
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        h_rms_layer[bin_idx-1] = std::vector<std::vector<std::shared_ptr<TH1D>>> (lambda_values.num+1);
        for (int lambda_idx=0; lambda_idx<=lambda_values.num; ++lambda_idx)
            h_rms_layer[bin_idx-1][lambda_idx] = std::vector<std::shared_ptr<TH1D>> (DAMPE_bgo_nLayers);
    }

    return h_rms_layer;
}