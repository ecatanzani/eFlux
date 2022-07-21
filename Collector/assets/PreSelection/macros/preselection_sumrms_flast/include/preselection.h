#ifndef PRESELECTION_H
#define PRESELECTION_H

extern void preselection(
        const char* input_list,
        const char* output_file,
        const bool verbose,
        const bool mc_flag,
        const unsigned int threads,
        const char* energy_config_file,
        const double input_spectral_index);

#endif