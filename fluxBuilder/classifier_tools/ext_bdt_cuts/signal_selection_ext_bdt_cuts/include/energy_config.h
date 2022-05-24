#ifndef ENERGY_CONFIG_H
#define ENERGY_CONFIG_H

#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

class energy_config
{
public:
    energy_config(const char* energy_config_file);
    ~energy_config(){};
    std::vector<double> GetEnergyBinning();
    const double GetMinEvtEnergy();
    const double GetMaxEvtEnergy();
    void PrintActiveFilters();

private:
    std::string parse_config_file(const char* config_file_path);
    void get_config_info(const std::string parsed_config);

    std::size_t n_bins;
    double min_event_energy {-999};
    double max_event_energy {-999};
    std::vector<double> energy_binning;
};

#endif