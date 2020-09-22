#include "BGO_energy_cuts.h"

void evaluateEnergyRatio(
    const std::vector<double> fracLayer,
    TH1D &h_layer_max_energy_ratio,
    std::vector<TH1D> &h_layer_energy_ratio,
    const double energy_w)
{
    for(auto it = fracLayer.begin(); it != fracLayer.end(); ++it)
        h_layer_energy_ratio[std::distance(fracLayer.begin(), it)].Fill(*it, energy_w);
    auto it_max = std::max_element(fracLayer.begin(), fracLayer.end());
    h_layer_max_energy_ratio.Fill(fracLayer[std::distance(fracLayer.begin(), it_max)], energy_w);
}