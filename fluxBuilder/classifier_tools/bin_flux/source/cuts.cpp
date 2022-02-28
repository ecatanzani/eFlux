#include "cuts.h"

cuts::cuts(
    const size_t n_cuts,
    const double bdt_cut_le,
    const double bdt_cut_he,
    const double xtrl_cut_le,
    const double xtrl_cut_he) {

        bdt_cuts.resize(n_cuts);
        xtrl_cuts.resize(n_cuts);
        bdt_xtrl_cuts.resize(pow(n_cuts, 2));

        for (size_t idx = 0; idx < n_cuts; ++idx) {
            double bdt_cut = bdt_cut_le + (bdt_cut_he - bdt_cut_le) * idx / (n_cuts - 1);
            double xtrl_cut = xtrl_cut_le + (xtrl_cut_he - xtrl_cut_le) * idx / (n_cuts - 1);
            bdt_cuts[idx] = std::make_tuple(bdt_cut, 0);
            xtrl_cuts[idx] = std::make_tuple(xtrl_cut, 0);
        }

        size_t idx = 0;
        for (size_t b_idx = 0; b_idx < n_cuts; ++b_idx) {
            for (size_t x_idx = 0; x_idx < n_cuts; ++x_idx) {
                bdt_xtrl_cuts[idx] = std::make_tuple(
                    std::get<0>(bdt_cuts[b_idx]),
                    std::get<0>(xtrl_cuts[x_idx]),
                    0);
                    ++idx;
            }
        }
    }

std::vector<std::tuple<double, unsigned int>> cuts::GetBDT() {
    return bdt_cuts;
}

std::vector<std::tuple<double, unsigned int>> cuts::GetXTRL() {
    return xtrl_cuts;
}

std::vector<std::tuple<double, double, unsigned int>> cuts::GetBDT_XTRL() {
    return bdt_xtrl_cuts;
}