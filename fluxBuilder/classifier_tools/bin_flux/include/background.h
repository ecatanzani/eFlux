#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <vector>
#include <string>
#include <tuple>

extern const std::vector<std::tuple<double, double>> estimate_background(
    std::string background_linear_fit_file,
    const std::vector<std::tuple<double, unsigned int>> bdt_cuts,
    const unsigned int energy_bin,
    const bool verbose,
    const unsigned int threads);

#endif