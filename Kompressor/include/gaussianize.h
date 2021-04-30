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

extern void studyGaussianizedTMVAvars(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    std::shared_ptr<lambda_config> _lambda_config,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE,
    const unsigned int threads,
    const bool _mc);

extern void fitGaussianizedTMVAvars(
    const std::string input_file,
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
    const rms_lambdas lambda_values,
    const int DAMPE_bgo_nLayers);

extern std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetFracLayerHistos(
    const int energy_nbins, 
    const energylastfraction_lambdas lambda_values,
    const int DAMPE_bgo_nLayers);

extern std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetRMSLayerHistosPointers(
    const int energy_nbins, 
    const rms_lambdas lambda_values,
    const int DAMPE_bgo_nLayers);

extern std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> GetAutoRMSLayerHistos(
    const int energy_nbins, 
    const rms_lambdas lambda_values,
    const int DAMPE_bgo_nLayers,
    const std::vector<std::vector<std::vector<double>>> &rms_boundaries);

extern void ComputeGoodness(std::vector<std::shared_ptr<TH1D>> h_rms);

extern void extract_lamda_info(
    const std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> input_histo,
    std::vector<double> &goodness,
    std::vector<double> &best_lambda,
    std::vector<int> &best_lambda_idx,
    const double lambda_start,
    const double lambda_step,
    const int bin_idx,
    const int lambda_idx,
    const int layer);

#endif