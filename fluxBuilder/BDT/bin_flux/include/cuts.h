#ifndef CUTS_H
#define CUTS_H

#include <cmath>
#include <tuple>
#include <vector>

class cuts {
    public:
        cuts(
            const size_t n_cuts = 1000,
            const double bdt_cut_le = -0.6,
            const double bdt_cut_he = 0.6,
            const double xtrl_cut_le = 0,
            const double xtrl_cut_he = 100);

        ~cuts() {};

        std::vector<std::tuple<double, unsigned int>> GetBDT();
        std::vector<std::tuple<double, unsigned int>> GetXTRL();
        std::vector<std::tuple<double, double, unsigned int>> GetBDT_XTRL();

    private:
        std::vector<std::tuple<double, unsigned int>> bdt_cuts;
        std::vector<std::tuple<double, unsigned int>> xtrl_cuts;
        std::vector<std::tuple<double, double, unsigned int>> bdt_xtrl_cuts;

};

#endif