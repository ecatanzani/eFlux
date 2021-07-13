#ifndef ENERGY_CONFIG_H
#define ENERGY_CONFIG_H

#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

struct tmva_set_energy_filter
{
    double min_event_energy = -999;
    double max_event_energy = -999;
};
class energy_config
{
public:
    energy_config(const std::string working_dir);
    ~energy_config(){};
    std::vector<float> GetEnergyBinning();
    void PrintActiveFilters();
    const double GetSetMinEvtEnergy();
    const double GetSetMaxEvtEnergy();

private:
    std::string parse_config_file(
        const std::string wd,
        const std::string config_file);
    void get_config_info(const std::string parsed_config);

    std::size_t n_bins;
    double min_event_energy = -999;
    double max_event_energy = -999;
    std::vector<float> energy_binning;
    tmva_set_energy_filter set_filter;
};

#endif