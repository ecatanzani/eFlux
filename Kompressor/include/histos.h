#ifndef HISTOS_H
#define HISTOS_H

#include <memory>
#include <vector>

#include "DAMPE_geo_structure.h"
#include "lambda_config.h"

#include "TH1D.h"
#include "TFile.h"

extern std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetRMSLayerHistos(
    const int energy_nbins, 
    const rms_lambdas lambda_values,
    const int DAMPE_bgo_nLayers);

extern std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetAutoRMSLayerHistos(
    const int energy_nbins, 
    const rms_lambdas lambda_values,
    const int DAMPE_bgo_nLayers,
    const std::vector<std::vector<std::vector<double>>> &rms_boundaries);

extern std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetFracLayerHistos(
    const int energy_nbins, 
    const energylastfraction_lambdas lambda_values,
    const int DAMPE_bgo_nLayers);

extern std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetHistos(
    TFile* infile,
    const int energy_nbins, 
    const int lambda_values_num,
    const double lambda_values_start,
    const double lambda_values_step,
    const int DAMPE_bgo_nLayers,
    const std::string abs_path);

#endif