#include "mc_ancillary.h"

void compute_proton_background(
    const double sumRms,
    const std::vector<double> fracLayer,
    const cuts_conf data_cuts,
    const double energy,
    TH1D &h_background_under_xtrl_cut,
    TH1D &h_background_over_xtrl_cut)
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

    if (xtrl < data_cuts.xtrl)
        h_background_under_xtrl_cut.Fill(energy * _GeV);
    if (xtrl > 20 && xtrl < 100)
        h_background_over_xtrl_cut.Fill(energy * _GeV);

}