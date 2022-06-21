#ifndef LATERAL_H
#define LATERAL_H

extern void buildLateralDistributions(
        const char* input_list,
        const char* output_file,
        const bool verbose,
        const bool mc_flag,
        const unsigned int threads,
        const std::string energy_config_file,
        const double input_spectral_index);

#endif