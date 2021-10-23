#ifndef HISTOS_H
#define HISTOS_H

#include "TH1D.h"
#include "TH2D.h"

#include <memory>
#include <vector>
#include <string>

#include "energy_config.h"

class histos {
    public:
        histos(std::shared_ptr<energy_config> econfig);
        ~histos() {};

        void Write(const std::string output_wd, const bool verbose);
    
        std::shared_ptr<TH1D> h_energy_fraction;
        std::shared_ptr<TH1D> h_energy_fraction_no_trigger;

        std::shared_ptr<TH1D> h_energy_fraction_large_angles_56;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_60;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_61;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_62;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_63;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_64;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_65;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_70;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_80;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_85;

        std::shared_ptr<TH2D> h_bar_energy;
        std::shared_ptr<TH1D> h_bars_last_layer_10MeV;
        std::shared_ptr<TH2D> h_bars_last_layer_10MeV_2D;

        std::shared_ptr<TH1D> h_maxrms;
        std::shared_ptr<TH2D> h_maxrms_2D;
        std::shared_ptr<TH1D> h_maxrms_no_trigger;
        std::shared_ptr<TH2D> h_maxrms_2D_no_trigger;

        std::shared_ptr<TH1D> h_bgoshower_top_X;
        std::shared_ptr<TH1D> h_bgoshower_bottom_X;
        std::shared_ptr<TH1D> h_bgoshower_top_Y;
        std::shared_ptr<TH1D> h_bgoshower_bottom_Y;

        std::shared_ptr<TH1D> h_energy_fraction_sh_axis_contained;
        std::shared_ptr<TH1D> h_energy_fraction_sh_axis_not_contained;

        private:
            int energy_nbins {0};
            std::vector<float> energy_binning;
            std::vector<float> bars_energy_bins;
            std::vector<float> number_of_bars_last_layer;
            std::vector<float> max_rms_bins;

};

#endif