#include "gaussianize.h"

#include <map>
#include <cmath>
#include <numeric>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>

#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"

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
    /*
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
    */
   
    #if 0
    // Create output log file
    std::ofstream _log("lambdafit_boxcox.log");
    std::vector<std::vector<std::shared_ptr<TH1D>>> best_rmshist (energy_nbins, std::vector<std::shared_ptr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<std::vector<double>> rms_lambda_corrections (energy_nbins);

    // Compute the goodness
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        std::vector<double> goodness_rms_layer (DAMPE_bgo_nLayers, 999);
        std::vector<double> goodness_frac_layer (DAMPE_bgo_nLayers, 999);

        std::vector<double> rms_best_lambda (DAMPE_bgo_nLayers, lambda_values.start);
        std::vector<double> fraclayer_best_lambda (DAMPE_bgo_nLayers, lambda_values.start);
        std::vector<int> rms_best_lambda_idx (DAMPE_bgo_nLayers, 0);
        std::vector<int> fraclayer_best_lambda_idx (DAMPE_bgo_nLayers, 0);

        int rms_statistics = h_rmsLayer_gauss[bin_idx-1][0][0]->GetEntries();
        //int fraclast_statistics = h_rmsLayer_gauss[bin_idx-1][0][0]->GetEntries();

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
                    
                    if (dof)
                    {
                        double tmp_goodness = chi2/dof;
                        if (goodness_rms_layer[ly]==999)
                        {
                            goodness_rms_layer[ly] = tmp_goodness;
                            best_lambda[ly] = lambda_values.start + lambda_values.step*lambda_idx;
                            best_lambda_idx[ly] = lambda_idx;
                        }
                        else
                        {
                            if (abs(tmp_goodness-1) < abs(goodness_rms_layer[ly]-1))
                            {
                                goodness_rms_layer[ly] = tmp_goodness;
                                best_lambda[ly] = lambda_values.start + lambda_values.step*lambda_idx;
                                best_lambda_idx[ly] = lambda_idx;
                            }
                        } 
                    }
                }
                rms_lambda_corrections[bin_idx-1] = best_lambda;
            }
            
            std::cout << "\n[ENERGY BIN " << bin_idx << "]\t - BGO Layer " << ly+1 << "\t - best lambda value: " << best_lambda[ly] << "\t - goodness: " << 
            goodness_rms_layer[ly] << "\t Stat (#events) [" << statistics << "]" << "\t best lambda idx: " << best_lambda_idx[ly];
            _log << "\n[ENERGY BIN " << bin_idx << "]\t - BGO Layer " << ly+1 << "\t - best lambda value: " << best_lambda[ly] << "\t - goodness: " << goodness_rms_layer[ly] << "\t Stat (#events) [" << statistics << "]" << "\t best lambda idx: " << best_lambda_idx[ly];
            best_rmshist[bin_idx-1][ly] = h_rmsLayer_gauss[bin_idx-1][best_lambda_idx[ly]][ly];
        }
    }

    // Write output file
    TFile* output_file = TFile::Open(outputPath.c_str(), "RECREATE");
    if (output_file->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outputPath << "]\n\n";
        exit(100);
    }

    // Write histos and graphs
    // Building lambda vs energy plots for each BGO layers
    std::vector<double> ebins (energy_nbins);
    std::vector<std::unique_ptr<TGraph>> lambda_rms_gr (DAMPE_bgo_nLayers);
    std::iota (std::begin(ebins), std::end(ebins), 1);

    auto get_lambda_per_layer = [&lambda_values](const std::vector<std::vector<double>> in_lambda, const int layer) -> const std::vector<double>
    {
        std::vector<double> out_lambda (in_lambda.size());
        if (layer<DAMPE_bgo_nLayers)
            for (unsigned int idx=0; idx<in_lambda.size(); ++idx)
                out_lambda[idx] = in_lambda[idx][layer] != lambda_values.start ? in_lambda[idx][layer] : 0;
        else
            for (auto&& elm:out_lambda)
                elm = 0;
        return out_lambda;
    };

    for (int ly=0; ly<DAMPE_bgo_nLayers; ++ly)
    {
        lambda_rms_gr[ly] = std::make_unique<TGraph> (energy_nbins, &ebins[0], &(get_lambda_per_layer(rms_lambda_corrections, ly))[0]);
        lambda_rms_gr[ly]->SetName((std::string("rms_lambda_energyevolution_layer_") + std::to_string(ly+1)).c_str());
        lambda_rms_gr[ly]->SetTitle((std::string("RMS #lambda - layer ") + std::to_string(ly+1)).c_str());
        lambda_rms_gr[ly]->GetXaxis()->SetTitle("energy layer");
        lambda_rms_gr[ly]->GetYaxis()->SetTitle("#lambda");
        lambda_rms_gr[ly]->Fit("pol1", "Q");
        lambda_rms_gr[ly]->Write();
    }

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        output_file->mkdir((std::string("energybin_") + std::to_string(bin_idx)).c_str());
        output_file->cd((std::string("energybin_") + std::to_string(bin_idx)).c_str());
        for (int lambda_idx=0; lambda_idx<=lambda_values.num; ++lambda_idx)
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->Write();
    }

    // Create final canvases
    std::vector<std::shared_ptr<TCanvas>> c_bestfit (energy_nbins);

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        c_bestfit[bin_idx-1] = std::make_shared<TCanvas> ((std::string("energybin_") + std::to_string(bin_idx)).c_str(), (std::string("energybin_") + std::to_string(bin_idx)).c_str());
        c_bestfit[bin_idx-1]->Divide(7,2);
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            c_bestfit[bin_idx-1]->cd(ly+1);
            best_rmshist[bin_idx-1][ly]->Draw();
        }
        output_file->cd((std::string("energybin_") + std::to_string(bin_idx)).c_str());
        c_bestfit[bin_idx-1]->cd(0);
        c_bestfit[bin_idx-1]->Write();
    }

    output_file->Close();

    // Writing TTree lambda correction file
    std::string tree_vars_file = outputPath.substr(0, outputPath.find_last_of("/")) + std::string("/lambda_corrections.root");
    TFile* corrections_file = TFile::Open(tree_vars_file.c_str(), "RECREATE");
    if (corrections_file->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << tree_vars_file << "]\n\n";
        exit(100);
    }

    TTree corrections_tree("corrections_tree", "Lambda corrections");
    
    int energy_bin_idx = 0;
    std::vector<double> rms_lambda (DAMPE_bgo_nLayers);

    corrections_tree.Branch("energy_bin_idx", &energy_bin_idx, "energy_bin_idx/I");
    corrections_tree.Branch("rms_lambda", &rms_lambda);

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        energy_bin_idx = bin_idx;
        rms_lambda = rms_lambda_corrections[bin_idx-1];
        corrections_tree.Fill();
    }

    corrections_tree.Write();
    corrections_file->Close();

    if (_VERBOSE)
    {
        std::cout << "\n\nOutput file has been written... [" << outputPath << "]\n";
        std::cout << "Output corrections file has been written... [" << tree_vars_file << "]\n";
    }

    #endif
}

std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetRMSLayerHistosPointers(
    const int energy_nbins, 
    const rms_lambdas lambda_values,
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