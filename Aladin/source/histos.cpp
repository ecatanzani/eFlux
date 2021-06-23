#include "histos.h"

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
    int nbins = 1000;
    double hmin = 0;
    double hmax = 200;
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
                std::string str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
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
    double hmax = 100;
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
                std::string str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
                std::string h_name = std::string("h_fraclayer_energybin_") + std::to_string(bin_idx) + std::string("_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
                h_frac_layer[bin_idx-1][lambda_idx][ly] = std::make_shared<TH1D>(h_name.c_str(), h_name.c_str(), nbins, hmin, hmax);
            }
        }
    }
    return h_frac_layer;
}