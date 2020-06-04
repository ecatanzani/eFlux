#ifndef ENERGY_CUTS_H
#define ENERGY_CUTS_H

#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

struct energy_cuts
{
    double cut_min_energy;
    double cut_max_energy;
    unsigned int e_bins;
};

extern void load_energy_struct(energy_cuts &e_cuts, const std::string wd);

#endif