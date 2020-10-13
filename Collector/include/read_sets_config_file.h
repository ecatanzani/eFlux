#ifndef READ_SETS_CONFIG_FILE_H
#define READ_SETS_CONFIG_FILE_H

#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

struct data_set_conf
{
    std::vector<std::string> sets_name;
    std::vector<double> sets_eMin;
    std::vector<double> sets_eMax;
    std::vector<int> sets_powerlaw_idx;
};

extern void load_input_dsets_config(data_set_conf &input_sets, const std::string wd);
extern int getInputPowerLawIndex(const double lEnergy, const double hEnergy, const data_set_conf input_sets);

#endif