#ifndef GAUSSIANIZE_H
#define GAUSSIANIZE_H

#include <string>
#include <memory>
#include <vector>

#include "config.h"
#include "energy_config.h"
#include "lambda_config.h"
#include "DAMPE_geo_structure.h"

#include "TH1D.h"
#include "TChain.h"
#include <ROOT/RDataFrame.hxx>

extern void gaussianizeTMVAvars(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    std::shared_ptr<lambda_config> _lambda_config,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE,
    const unsigned int threads,
    const bool _mc);

extern void showGaussianizedTMVAvars(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    std::shared_ptr<lambda_config> _lambda_config,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE,
    const unsigned int threads,
    const bool _mc);

extern std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetRMSLayerHistos(
    const int energy_nbins, 
    const lambdas lambda_values,
    const int DAMPE_bgo_nLayers);

extern std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetAutoRMSLayerHistos(
    const int energy_nbins, 
    const lambdas lambda_values,
    const int DAMPE_bgo_nLayers,
    const std::vector<std::vector<std::vector<double>>> &rms_boundaries);

extern void ComputeGoodness(std::vector<std::shared_ptr<TH1D>> h_rms);

#endif