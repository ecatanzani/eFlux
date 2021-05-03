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
    // Extract the energy binning
    auto energy_binning = _energy_config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;
    // Extract lambda config
    auto rms_lambda_values = _lambda_config->GetRMSLambdaStruct();
    auto elf_lambda_values = _lambda_config->GetELFLambdaStruct();
    if (_VERBOSE) _lambda_config->PrintLambdaSettings();
    
    // Read input file
    TFile *infile = TFile::Open(input_file.c_str(), "READ");
    if (infile->IsZombie())
    {
        std::cerr << "\n\nError opening input file [" << input_file << "]\n\n";
        exit(100);
    }

    // Get histo structure
    auto h_rmsLayer_gauss = GetHistos(
        infile, 
        energy_nbins, 
        rms_lambda_values.num, 
        rms_lambda_values.start,
        rms_lambda_values.step, 
        DAMPE_bgo_nLayers,
        "RMS");
    auto h_fracLayer_gauss = GetHistos(
        infile, 
        energy_nbins, 
        elf_lambda_values.num, 
        elf_lambda_values.start,
        elf_lambda_values.step, 
        DAMPE_bgo_nLayers,
        "ELF");
    
    infile->Close();
    
    // Create output log file
    std::ofstream _log("lambdafit_boxcox.log");
    std::vector<std::vector<std::shared_ptr<TH1D>>> best_rmshist (energy_nbins, std::vector<std::shared_ptr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<std::vector<std::shared_ptr<TH1D>>> best_fraclasthist (energy_nbins, std::vector<std::shared_ptr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<std::vector<double>> rms_lambda_corrections (energy_nbins);
    std::vector<std::vector<double>> fraclast_lambda_corrections (energy_nbins);

    // Compute the goodness
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        std::vector<double> goodness_rms_layer (DAMPE_bgo_nLayers, 999);
        std::vector<double> goodness_frac_layer (DAMPE_bgo_nLayers, 999);

        std::vector<double> rms_best_lambda (DAMPE_bgo_nLayers, rms_lambda_values.start);
        std::vector<double> fraclayer_best_lambda (DAMPE_bgo_nLayers, elf_lambda_values.start);
        std::vector<int> rms_best_lambda_idx (DAMPE_bgo_nLayers, 0);
        std::vector<int> fraclayer_best_lambda_idx (DAMPE_bgo_nLayers, 0);

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
                    elf_lambda_values.start,
                    elf_lambda_values.step,
                    bin_idx, 
                    lambda_idx, 
                    ly);
        
        rms_lambda_corrections[bin_idx-1] = rms_best_lambda;
        fraclast_lambda_corrections[bin_idx-1] = fraclayer_best_lambda;

        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            best_rmshist[bin_idx-1][ly] = h_rmsLayer_gauss[bin_idx-1][rms_best_lambda_idx[ly]][ly];
            best_fraclasthist[bin_idx-1][ly] = h_fracLayer_gauss[bin_idx-1][fraclayer_best_lambda_idx[ly]][ly];
        }

        if (_VERBOSE)
        {
            std::cout << "\n\nEnergy bin [" << bin_idx << "]\n";
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                std::cout << "\n(RMS) - BGO Layer " << ly+1 << "\t - best lambda value: " << rms_best_lambda[ly] << "\t - goodness: " << 
            goodness_rms_layer[ly] << "\t Stat (#events) [" << rms_statistics << "]" << "\t best lambda idx: " << rms_best_lambda_idx[ly];
                _log << "\n(RMS) - BGO Layer " << ly+1 << "\t - best lambda value: " << rms_best_lambda[ly] << "\t - goodness: " << 
            goodness_rms_layer[ly] << "\t Stat (#events) [" << rms_statistics << "]" << "\t best lambda idx: " << rms_best_lambda_idx[ly];
            }
            std::cout << std::endl;
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                std::cout << "\n(EnergyLastLayerFraction) - BGO Layer " << ly+1 << "\t - best lambda value: " << fraclayer_best_lambda[ly] << "\t - goodness: " << 
            goodness_frac_layer[ly] << "\t Stat (#events) [" << fraclast_statistics << "]" << "\t best lambda idx: " << fraclayer_best_lambda_idx[ly];
                _log << "\n(EnergyLastLayerFraction) - BGO Layer " << ly+1 << "\t - best lambda value: " << fraclayer_best_lambda[ly] << "\t - goodness: " << 
            goodness_frac_layer[ly] << "\t Stat (#events) [" << fraclast_statistics << "]" << "\t best lambda idx: " << fraclayer_best_lambda_idx[ly];
            }
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
    std::vector<std::unique_ptr<TGraph>> lambda_fraclast_gr (DAMPE_bgo_nLayers);
    std::iota (std::begin(ebins), std::end(ebins), 1);

    auto get_lambda_per_layer = [](const std::vector<std::vector<double>> in_lambda, const int layer, const double lambda_start) -> const std::vector<double>
    {
        std::vector<double> out_lambda (in_lambda.size());
        if (layer<DAMPE_bgo_nLayers)
            for (unsigned int idx=0; idx<in_lambda.size(); ++idx)
                out_lambda[idx] = in_lambda[idx][layer] != lambda_start ? in_lambda[idx][layer] : 0;
        else
            for (auto&& elm:out_lambda)
                elm = 0;
        return out_lambda;
    };

    for (int ly=0; ly<DAMPE_bgo_nLayers; ++ly)
    {
        // Build RMS graph
        lambda_rms_gr[ly] = std::make_unique<TGraph> (energy_nbins, &ebins[0], &(get_lambda_per_layer(rms_lambda_corrections, ly, rms_lambda_values.start))[0]);
        lambda_rms_gr[ly]->SetName((std::string("rms_lambda_energyevolution_layer_") + std::to_string(ly+1)).c_str());
        lambda_rms_gr[ly]->SetTitle((std::string("RMS #lambda - layer ") + std::to_string(ly+1)).c_str());
        lambda_rms_gr[ly]->GetXaxis()->SetTitle("energy layer");
        lambda_rms_gr[ly]->GetYaxis()->SetTitle("#lambda");
        lambda_rms_gr[ly]->Fit("pol1", "Q");
        lambda_rms_gr[ly]->Write();
        // Build energy fraction graph
        lambda_fraclast_gr[ly] = std::make_unique<TGraph> (energy_nbins, &ebins[0], &(get_lambda_per_layer(fraclast_lambda_corrections, ly, elf_lambda_values.start))[0]);
        lambda_fraclast_gr[ly]->SetName((std::string("fraclast_lambda_energyevolution_layer_") + std::to_string(ly+1)).c_str());
        lambda_fraclast_gr[ly]->SetTitle((std::string("Energy Fraction #lambda - layer ") + std::to_string(ly+1)).c_str());
        lambda_fraclast_gr[ly]->GetXaxis()->SetTitle("energy layer");
        lambda_fraclast_gr[ly]->GetYaxis()->SetTitle("#lambda");
        lambda_fraclast_gr[ly]->Fit("pol1", "Q");
        lambda_fraclast_gr[ly]->Write();
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
    
    // Create final canvases
    std::vector<std::shared_ptr<TCanvas>> rms_c_bestfit (energy_nbins);
    std::vector<std::shared_ptr<TCanvas>> fraclast_c_bestfit (energy_nbins);

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        // RMS canvasis
        rms_c_bestfit[bin_idx-1] = std::make_shared<TCanvas> ((std::string("energybin_") + std::to_string(bin_idx)).c_str(), (std::string("energybin_") + std::to_string(bin_idx)).c_str());
        rms_c_bestfit[bin_idx-1]->Divide(7,2);
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            rms_c_bestfit[bin_idx-1]->cd(ly+1);
            best_rmshist[bin_idx-1][ly]->Draw();
        }
        output_file->cd((std::string("RMS/energybin_") + std::to_string(bin_idx)).c_str());
        rms_c_bestfit[bin_idx-1]->cd(0);
        rms_c_bestfit[bin_idx-1]->Write();

        // Energy fraction canvases
        fraclast_c_bestfit[bin_idx-1] = std::make_shared<TCanvas> ((std::string("energybin_") + std::to_string(bin_idx)).c_str(), (std::string("energybin_") + std::to_string(bin_idx)).c_str());
        fraclast_c_bestfit[bin_idx-1]->Divide(7,2);
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            fraclast_c_bestfit[bin_idx-1]->cd(ly+1);
            best_fraclasthist[bin_idx-1][ly]->Draw();
        }
        output_file->cd((std::string("ELF/energybin_") + std::to_string(bin_idx)).c_str());
        fraclast_c_bestfit[bin_idx-1]->cd(0);
        fraclast_c_bestfit[bin_idx-1]->Write();
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
    std::vector<double> energyfractionll_lambda (DAMPE_bgo_nLayers);

    corrections_tree.Branch("energy_bin_idx", &energy_bin_idx, "energy_bin_idx/I");
    corrections_tree.Branch("rms_lambda", &rms_lambda);
    corrections_tree.Branch("energyfractionll_lambda", &energyfractionll_lambda);

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        energy_bin_idx = bin_idx;
        rms_lambda = rms_lambda_corrections[bin_idx-1];
        energyfractionll_lambda = fraclast_lambda_corrections[bin_idx-1];
        corrections_tree.Fill();
    }

    corrections_tree.Write();
    corrections_file->Close();

    if (_VERBOSE)
    {
        std::cout << "\n\nOutput file has been written... [" << outputPath << "]\n";
        std::cout << "Output corrections file has been written... [" << tree_vars_file << "]\n";
    }
    
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

std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetHistos(
    TFile* infile,
    const int energy_nbins, 
    const int lambda_values_num,
    const double lambda_values_start,
    const double lambda_values_step,
    const int DAMPE_bgo_nLayers,
    const std::string abs_path)
{
    // Get automatically the histo name
    auto tmphname = !strcmp(abs_path.c_str(), "RMS") ? std::string("h_rms_energybin_") : std::string("h_fraclayer_energybin_");
    // Build & Fill the histo vector stucture
    std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> histos (energy_nbins);
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        auto lambda = lambda_values_start;
        histos[bin_idx-1] = std::vector<std::vector<std::shared_ptr<TH1D>>> (lambda_values_num+1);
        for (int lambda_idx=0; lambda_idx<=lambda_values_num; ++lambda_idx)
        {
            histos[bin_idx-1][lambda_idx] = std::vector<std::shared_ptr<TH1D>> (DAMPE_bgo_nLayers);
            if (lambda_idx) lambda += lambda_values_step;
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                std::string str_lambda = lambda<0 ? std::string("neg_") + std::to_string(abs(lambda)) : std::to_string(lambda);
                std::string h_name = abs_path + std::string("/energybin_") + std::to_string(bin_idx) + std::string("/") + tmphname + std::to_string(bin_idx) + std::string("_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
                histos[bin_idx-1][lambda_idx][ly] = std::shared_ptr<TH1D>(static_cast<TH1D*>(infile->Get((h_name).c_str())));
                histos[bin_idx-1][lambda_idx][ly]->SetDirectory(0);
            }
        }
    }
    return histos;
}
