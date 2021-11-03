#ifndef HISTOS_H
#define HISTOS_H

#include "TH1D.h"
#include "TH2D.h"

#include <memory>
#include <vector>
#include <string>

#include "energy_config.h"

#include "DmpEvtSimuPrimaries.h"

class histos {
    public:
        histos(std::shared_ptr<energy_config> econfig, const bool mc);
        ~histos() {};
        void SetWeight(std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries, const double evt_corr_energy_gev);
        const double GetWeight();

        void Write(const std::string output_wd, const bool verbose);
    
        std::shared_ptr<TH1D> h_energy_fraction;
        std::shared_ptr<TH2D> h_energy_fraction_2D;
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

        std::shared_ptr<TH1D> h_energy_fraction_after_bgofiducial;
        std::shared_ptr<TH2D> h_energy_fraction_2D_after_bgofiducial;
        std::shared_ptr<TH1D> h_energy_fraction_no_trigger_after_bgofiducial;

        std::shared_ptr<TH1D> h_energy_fraction_large_angles_56_after_bgofiducial;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_60_after_bgofiducial;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_61_after_bgofiducial;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_62_after_bgofiducial;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_63_after_bgofiducial;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_64_after_bgofiducial;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_65_after_bgofiducial;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_70_after_bgofiducial;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_80_after_bgofiducial;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_85_after_bgofiducial;

        std::shared_ptr<TH1D> h_energy_fraction_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_energy_fraction_2D_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_energy_fraction_no_trigger_bgofiducial_lastcut;

        std::shared_ptr<TH1D> h_energy_fraction_large_angles_56_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_60_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_61_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_62_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_63_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_64_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_65_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_70_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_80_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_energy_fraction_large_angles_85_bgofiducial_lastcut;

        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_0;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_1;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_2;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_3;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_4;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_5;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_6;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_7;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_8;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_9;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_10;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_11;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_12;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_13;

        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_0_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_1_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_2_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_3_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_4_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_5_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_6_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_7_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_8_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_9_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_10_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_11_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_12_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_13_bgoshower_in;

        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_0_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_1_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_2_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_3_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_4_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_5_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_6_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_7_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_8_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_9_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_10_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_11_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_12_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_13_bgoshower_out;

        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial;

        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial_bgoshower_in;

        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial_bgoshower_out;

        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut;

        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut_bgoshower_in;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut_bgoshower_in;

        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut_bgoshower_out;
        std::shared_ptr<TH2D> h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut_bgoshower_out;

        std::shared_ptr<TH2D> h_bar_energy;
        std::shared_ptr<TH2D> h_bar_energy_2D;
        std::shared_ptr<TH1D> h_bars_last_layer_10MeV;
        std::shared_ptr<TH2D> h_bars_last_layer_10MeV_2D;

        std::shared_ptr<TH2D> h_bar_energy_after_bgofiducial;
        std::shared_ptr<TH2D> h_bar_energy_2D_after_bgofiducial;
        std::shared_ptr<TH1D> h_bars_last_layer_10MeV_after_bgofiducial;
        std::shared_ptr<TH2D> h_bars_last_layer_10MeV_2D_after_bgofiducial;

        std::shared_ptr<TH2D> h_bar_energy_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_bar_energy_2D_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_bars_last_layer_10MeV_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_bars_last_layer_10MeV_2D_bgofiducial_lastcut;

        std::shared_ptr<TH2D> h_bar_energy_lateral_showering_lastcut;
        std::shared_ptr<TH2D> h_bar_energy_2D_lateral_showering_lastcut;
        std::shared_ptr<TH1D> h_bars_last_layer_10MeV_lateral_showering_lastcut;
        std::shared_ptr<TH2D> h_bars_last_layer_10MeV_2D_lateral_showering_lastcut;

        std::shared_ptr<TH1D> h_maxrms;
        std::shared_ptr<TH2D> h_maxrms_2D;
        std::shared_ptr<TH1D> h_maxrms_no_trigger;
        std::shared_ptr<TH2D> h_maxrms_2D_no_trigger;

        std::shared_ptr<TH1D> h_maxrms_after_bgofiducial;
        std::shared_ptr<TH2D> h_maxrms_2D_after_bgofiducial;
        std::shared_ptr<TH1D> h_maxrms_no_trigger_after_bgofiducial;
        std::shared_ptr<TH2D> h_maxrms_2D_no_trigger_after_bgofiducial;

        std::shared_ptr<TH1D> h_maxrms_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_maxrms_2D_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_maxrms_no_trigger_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_maxrms_2D_no_trigger_bgofiducial_lastcut;

        std::shared_ptr<TH1D> h_maxrms_lateral_showering_lastcut;
        std::shared_ptr<TH2D> h_maxrms_2D_lateral_showering_lastcut;

        std::shared_ptr<TH1D> h_bgoshower_top_X;
        std::shared_ptr<TH1D> h_bgoshower_bottom_X;
        std::shared_ptr<TH1D> h_bgoshower_top_Y;
        std::shared_ptr<TH1D> h_bgoshower_bottom_Y;

        std::shared_ptr<TH1D> h_bgoshower_top_X_after_bgofiducial;
        std::shared_ptr<TH1D> h_bgoshower_bottom_X_after_bgofiducial;
        std::shared_ptr<TH1D> h_bgoshower_top_Y_after_bgofiducial;
        std::shared_ptr<TH1D> h_bgoshower_bottom_Y_after_bgofiducial;

        std::shared_ptr<TH1D> h_bgoshower_top_X_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_bgoshower_bottom_X_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_bgoshower_top_Y_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_bgoshower_bottom_Y_bgofiducial_lastcut;

        std::shared_ptr<TH2D> h_bgoshower_top_X_simu_reco_energy_diff;
        std::shared_ptr<TH2D> h_bgoshower_bottom_X_simu_reco_energy_diff;
        std::shared_ptr<TH2D> h_bgoshower_top_Y_simu_reco_energy_diff;
        std::shared_ptr<TH2D> h_bgoshower_bottom_Y_simu_reco_energy_diff;

        std::shared_ptr<TH2D> h_bgoshower_top_X_simu_reco_energy_diff_after_bgofiducial;
        std::shared_ptr<TH2D> h_bgoshower_bottom_X_simu_reco_energy_diff_after_bgofiducial;
        std::shared_ptr<TH2D> h_bgoshower_top_Y_simu_reco_energy_diff_after_bgofiducial;
        std::shared_ptr<TH2D> h_bgoshower_bottom_Y_simu_reco_energy_diff_after_bgofiducial;

        std::shared_ptr<TH2D> h_bgoshower_top_X_simu_reco_energy_diff_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_bgoshower_bottom_X_simu_reco_energy_diff_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_bgoshower_top_Y_simu_reco_energy_diff_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_bgoshower_bottom_Y_simu_reco_energy_diff_bgofiducial_lastcut;

        std::shared_ptr<TH1D> h_diff_bgo_simu_particle_direction_X;
        std::shared_ptr<TH1D> h_diff_bgo_simu_particle_direction_Y;
        std::shared_ptr<TH1D> h_diff_bgo_simu_extr_top_position_X;
        std::shared_ptr<TH1D> h_diff_bgo_simu_extr_top_position_Y;

        std::shared_ptr<TH1D> h_diff_bgo_simu_particle_direction_X_after_bgofiducial;
        std::shared_ptr<TH1D> h_diff_bgo_simu_particle_direction_Y_after_bgofiducial;
        std::shared_ptr<TH1D> h_diff_bgo_simu_extr_top_position_X_after_bgofiducial;
        std::shared_ptr<TH1D> h_diff_bgo_simu_extr_top_position_Y_after_bgofiducial;

        std::shared_ptr<TH1D> h_diff_bgo_simu_particle_direction_X_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_diff_bgo_simu_particle_direction_Y_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_diff_bgo_simu_extr_top_position_X_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_diff_bgo_simu_extr_top_position_Y_bgofiducial_lastcut;

        std::shared_ptr<TH1D> h_energy_fraction_sh_axis_contained;
        std::shared_ptr<TH1D> h_energy_fraction_sh_axis_not_contained;

        std::shared_ptr<TH1D> h_energy_fraction_sh_axis_contained_after_bgofiducial;
        std::shared_ptr<TH1D> h_energy_fraction_sh_axis_not_contained_after_bgofiducial;

        std::shared_ptr<TH1D> h_energy_fraction_sh_axis_contained_bgofiducial_lastcut;
        std::shared_ptr<TH1D> h_energy_fraction_sh_axis_not_contained_bgofiducial_lastcut;

        std::shared_ptr<TH1D> h_PSD_STK_X_match_energy_int;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_energy_int;
        std::shared_ptr<TH1D> h_PSD_STK_X_match_100_250;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_100_250;
        std::shared_ptr<TH1D> h_PSD_STK_X_match_250_500;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_250_500;
        std::shared_ptr<TH1D> h_PSD_STK_X_match_500_1000;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_500_1000;
        std::shared_ptr<TH1D> h_PSD_STK_X_match_1000_5000;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_1000_5000;
        std::shared_ptr<TH1D> h_PSD_STK_X_match_5000;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_5000;

        std::shared_ptr<TH1D> h_PSD_STK_X_match_energy_int_psd_fiducial;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_energy_int_psd_fiducial;
        std::shared_ptr<TH1D> h_PSD_STK_X_match_100_250_psd_fiducial;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_100_250_psd_fiducial;
        std::shared_ptr<TH1D> h_PSD_STK_X_match_250_500_psd_fiducial;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_250_500_psd_fiducial;
        std::shared_ptr<TH1D> h_PSD_STK_X_match_500_1000_psd_fiducial;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_500_1000_psd_fiducial;
        std::shared_ptr<TH1D> h_PSD_STK_X_match_1000_5000_psd_fiducial;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_1000_5000_psd_fiducial;
        std::shared_ptr<TH1D> h_PSD_STK_X_match_5000_psd_fiducial;
        std::shared_ptr<TH1D> h_PSD_STK_Y_match_5000_psd_fiducial;

        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_20_100;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_100_250;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_250_500;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_500_1000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_1000_3000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_3000_5000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_5000;
        
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_bgofiducial;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_bgofiducial_20_100;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_bgofiducial_100_250;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_bgofiducial_250_500;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_bgofiducial_500_1000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_bgofiducial_1000_3000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_bgofiducial_3000_5000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_bgofiducial_5000;

        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_20_100_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_100_250_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_250_500_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_500_1000_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_1000_3000_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_3000_5000_bgofiducial_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_5000_bgofiducial_lastcut;

        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_remove_lateral_and_showering;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_20_100;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_100_250;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_250_500;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_500_1000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_1000_3000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_3000_5000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_5000;

        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_lateral_showering_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_20_100_lateral_showering_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_100_250_lateral_showering_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_250_500_lateral_showering_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_500_1000_lateral_showering_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_1000_3000_lateral_showering_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_3000_5000_lateral_showering_lastcut;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_5000_lateral_showering_lastcut;

        std::shared_ptr<TH1D> h_STK_X_clusters;
        std::shared_ptr<TH1D> h_STK_Y_clusters;
        std::shared_ptr<TH1D> h_STK_X_holes;
        std::shared_ptr<TH1D> h_STK_Y_holes;

        std::shared_ptr<TH2D> h_STK_X_clusters_vs_energy;
        std::shared_ptr<TH2D> h_STK_Y_clusters_vs_energy;
        std::shared_ptr<TH2D> h_STK_X_holes_vs_energy;
        std::shared_ptr<TH2D> h_STK_Y_holes_vs_energy;

        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_difference;
        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_X_difference;
        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_Y_difference;
        std::shared_ptr<TH1D> h_STK_BGO_track_angular_difference;

        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_difference_3_clusters;
        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_X_difference_3_clusters;
        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_Y_difference_3_clusters;
        std::shared_ptr<TH1D> h_STK_BGO_track_angular_difference_3_clusters;

        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_difference_4_clusters;
        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_X_difference_4_clusters;
        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_Y_difference_4_clusters;
        std::shared_ptr<TH1D> h_STK_BGO_track_angular_difference_4_clusters;

        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_difference_5_clusters;
        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_X_difference_5_clusters;
        std::shared_ptr<TH1D> h_STK_BGO_TOP_spatial_Y_difference_5_clusters;
        std::shared_ptr<TH1D> h_STK_BGO_track_angular_difference_5_clusters;

        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_track_selection;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_track_selection_20_100;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_track_selection_100_250;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_track_selection_250_500;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_track_selection_500_1000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_track_selection_1000_3000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_track_selection_3000_5000;
        std::shared_ptr<TH2D> h_BGOrec_sumRms_flast_after_track_selection_5000;

        std::shared_ptr<TH1D> h_STK_charge_X;
        std::shared_ptr<TH1D> h_STK_charge_Y;
        std::shared_ptr<TH1D> h_STK_charge;
        std::shared_ptr<TH2D> h_STK_charge_2D;

        std::shared_ptr<TH1D> h_STK_charge_X_3_clusters;
        std::shared_ptr<TH1D> h_STK_charge_Y_3_clusters;
        std::shared_ptr<TH1D> h_STK_charge_3_clusters;
        std::shared_ptr<TH2D> h_STK_charge_2D_3_clusters;

        std::shared_ptr<TH1D> h_STK_charge_X_4_clusters;
        std::shared_ptr<TH1D> h_STK_charge_Y_4_clusters;
        std::shared_ptr<TH1D> h_STK_charge_4_clusters;
        std::shared_ptr<TH2D> h_STK_charge_2D_4_clusters;

        std::shared_ptr<TH1D> h_STK_charge_X_5_clusters;
        std::shared_ptr<TH1D> h_STK_charge_Y_5_clusters;
        std::shared_ptr<TH1D> h_STK_charge_5_clusters;
        std::shared_ptr<TH2D> h_STK_charge_2D_5_clusters;

        std::shared_ptr<TH1D> h_PSD_charge_X;
        std::shared_ptr<TH1D> h_PSD_charge_Y;
        std::shared_ptr<TH1D> h_PSD_charge;
        std::shared_ptr<TH2D> h_PSD_charge_2D;
        std::shared_ptr<TH1D> h_PSD_sum_of_XY_charges;

        std::shared_ptr<TH2D> h_PSD_X_clusters;
        std::shared_ptr<TH2D> h_PSD_Y_clusters;

        std::shared_ptr<TH1D> h_STK_charge_X_nocut;
        std::shared_ptr<TH1D> h_STK_charge_Y_nocut;
        std::shared_ptr<TH1D> h_STK_charge_nocut;
        std::shared_ptr<TH2D> h_STK_charge_2D_nocut;

        std::shared_ptr<TH1D> h_PSD_charge_X_nocut;
        std::shared_ptr<TH1D> h_PSD_charge_Y_nocut;
        std::shared_ptr<TH1D> h_PSD_charge_nocut;
        std::shared_ptr<TH2D> h_PSD_charge_2D_nocut;
        std::shared_ptr<TH1D> h_PSD_sum_of_XY_charges_nocut;

        std::shared_ptr<TH1D> h_STK_charge_X_PSD_charge_cut;
        std::shared_ptr<TH1D> h_STK_charge_Y_PSD_charge_cut;
        std::shared_ptr<TH1D> h_STK_charge_PSD_charge_cut;
        std::shared_ptr<TH2D> h_STK_charge_2D_PSD_charge_cut;

        std::shared_ptr<TH1D> h_PSD_charge_X_STK_charge_cut;
        std::shared_ptr<TH1D> h_PSD_charge_Y_STK_charge_cut;
        std::shared_ptr<TH1D> h_PSD_charge_STK_charge_cut;
        std::shared_ptr<TH2D> h_PSD_charge_2D_STK_charge_cut;
        std::shared_ptr<TH1D> h_PSD_sum_of_XY_charges_STK_charge_cut;

        private:
            double weight {1};
            bool h_simu {false};
            int energy_nbins {0};

            std::vector<float> energy_binning;
            std::vector<float> energy_fraction_bins;
            std::vector<float> bars_energy_bins;
            std::vector<float> bars_2d_energy_bins;
            std::vector<float> number_of_bars_last_layer;
            std::vector<float> max_rms_bins;
            std::vector<float> sumRms_binning;
            std::vector<float> flast_binning;
            std::vector<float> stk_track_binning;
            std::vector<float> psd_clusters_binning;

};

#endif