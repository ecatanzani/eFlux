#ifndef BDT_CONFIG_H
#define BDT_CONFIG_H

#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

class bdt_config {
    public:
        bdt_config(const char* config_file);
        ~bdt_config(){};
        
        const double GetLowEnergyBDTCut();
        const double GetMidEnergyBDTCut();
        const double GetHighEnergyBDTCut();

    private:
        std::string parse_config_file(const char* config_file_path);
        void get_config_info(std::string parsed_config);
        
        double le_c_cut {0};
        double me_c_cut {0};
        double he_c_cut {0};
};

#endif