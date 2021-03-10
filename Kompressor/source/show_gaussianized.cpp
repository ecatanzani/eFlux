#include "gaussianize.h"

#include <map>
#include <cmath>
#include <vector>
#include <memory>

#include "DAMPE_geo_structure.h"

#include "TH1D.h"
#include <ROOT/RDataFrame.hxx>

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
    
    std::cout << "\n\n**** RDF statistics ****\n";
    std::cout << "***************************\n";
    std::cout << "\nPreselected events: " << *(_data_fr.Count());
    std::cout << "\n\n********************\n";

    auto lambda_values = _lambda_config->GetLambdaStruct();
    if (_VERBOSE)
        _lambda_config->PrintLambdaSettings();

    if (_VERBOSE)
        std::cout << "\nBuilding histos...";

    auto h_rmsLayer_gauss = GetRMSLayerHistos(
        energy_nbins, 
        lambda_values, 
        DAMPE_bgo_nLayers);

    if (_VERBOSE)
        std::cout << "\nAnalysis running...";

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {   
        auto bin_filter = [&bin_idx](int energy_bin) -> bool { return energy_bin == bin_idx; };
        _data_fr.Filter(bin_filter, {"energy_bin"})
                .Foreach([&h_rmsLayer_gauss, &bin_idx](const std::map<std::vector<double>, double> rmslayer_gauss, const double energy_w)
                {
                    auto lambda_idx = 0;
                    for (const auto &map_elm : rmslayer_gauss)
                    {
                        auto rms_gauss = map_elm.first;
                        //auto lambda = map_elm.second;
                        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                            h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->Fill(rms_gauss[ly], energy_w);
                        ++lambda_idx;
                    }
                }, {"rmsLayer_gauss", "simu_energy_w_corr"});
    }
    
    if (_VERBOSE)
        std::cout << "\n\nWriting output file... [" << outputPath << "]\n";
    
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
        for (int lambda_idx=0; lambda_idx<lambda_values.num; ++lambda_idx)
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                h_rmsLayer_gauss[bin_idx-1][lambda_idx][ly]->Write();
    }

    output_file->Close();

}

std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetRMSLayerHistos(
    const int energy_nbins, 
    const lambdas lambda_values,
    const int DAMPE_bgo_nLayers)
{
    int nbins = 1e+2;
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
                str_lambda.erase(str_lambda.find_last_not_of('0'), std::string::npos);
                std::string h_name = std::string("h_rms_energybin_") + std::to_string(bin_idx) + std::string("_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
                h_rms_layer[bin_idx-1][lambda_idx][ly] = std::make_shared<TH1D>(h_name.c_str(), h_name.c_str(), nbins, -1e-2, 1e-2);
            }
        }
    }
    return h_rms_layer;
}