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
    energy_config(const std::string working_dir);
    ~energy_config(){};
    const double GetMinEvtEnergy();
    const double GetMaxEvtEnergy();

private:
    std::string parse_config_file(const std::string wd, const std::string config_file);
    void get_config_info(const std::string parsed_config);
    
    double min_event_energy {-999};
    double max_event_energy {-999};
    const std::string config_file_name {"energy_config.conf"};
};

#endif