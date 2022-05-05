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
        
        /*
        const double GetLowEnergyBDTCut();
        const double GetMidEnergyBDTCut();
        const double GetHighEnergyBDTCut();
        */

        const double GetBDTCut_10_100();
        const double GetBDTCut_100_250();
        const double GetBDTCut_250_500();
        const double GetBDTCut_500_1000();
        const double GetBDTCut_1000_3000();
        const double GetBDTCut_3000();

    private:
        std::string parse_config_file(const char* config_file_path);
        void get_config_info(std::string parsed_config);
        
        /*
        double le_c_cut {0};
        double me_c_cut {0};
        double he_c_cut {0};
        */

        double cut_10_100 {0};
        double cut_100_250 {0};
        double cut_250_500 {0};
        double cut_500_1000 {0};
        double cut_1000_3000 {0};
        double cut_3000 {0};
};

#endif