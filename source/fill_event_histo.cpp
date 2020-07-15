#include "fill_event_histo.h"

void fill_XTRL_histo(
    const double sumRms,
    const std::vector<double> fracLayer,
    const cuts_conf data_cuts,
    const double energy,
    TH1D &h_xtrl_energy_int,
    TH2D &h_xtrl)
{
    unsigned int last_layer_idx = 0;
    double _GeV = 0.001;

    // Find the last layer with an energy release
    for (auto it = fracLayer.begin(); it != fracLayer.end(); ++it)
        if (*it)
        {
            auto index = std::distance(fracLayer.begin(), it);
            last_layer_idx = index;
        }

    // Build XTRL
    double xtrl = 0.1251e-8 * pow(sumRms, 4) * fracLayer[last_layer_idx];

    // Fill XTRL histos
    h_xtrl_energy_int.Fill(xtrl);
    h_xtrl.Fill(energy * _GeV, xtrl);
}