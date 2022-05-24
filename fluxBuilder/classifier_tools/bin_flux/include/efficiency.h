#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <tuple>
#include <string>
#include <vector>

#include "config.h"

extern const std::vector<std::tuple<double, double, double>> compute_efficiency(
    std::shared_ptr<config> bdt_config,
    const std::string learning_method,
    const std::vector<std::tuple<double, unsigned int>> bdt_cuts,
    const std::string mc_file_list,
    const char* mc_correction_function,
    const unsigned int energy_bin,
    std::vector<double>& energy_binning,
    const bool verbose,
    const unsigned int threads);

#endif