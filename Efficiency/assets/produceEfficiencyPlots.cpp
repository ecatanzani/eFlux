#include <iostream>
#include <string>

#include "TPDF.h"
#include "TPad.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAttMarker.h"
#include "TPaveLabel.h"
#include "TEfficiency.h"
#include "TLegendEntry.h"
#include "TGraphAsymmErrors.h"


void produceEfficiencyPlots(
    const char* data_file,
    const char* mc_file,
    const char* data_label="DATA",
    const char* mc_label="electron MC") {

        TFile* datafile = TFile::Open(data_file, "READ");
        if (datafile->IsZombie()) {
            std::cerr << "\n\nError opening input DATA ROOT file [" << data_file << "]\n\n";
            exit(100);
        }

        auto data_trigger_eff_het_over_let_xtrl_tight               = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_let_xtrl_tight"));
        auto data_trigger_eff_het_over_unb_xtrl_tight               = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_unb_xtrl_tight"));
        auto data_maxrms_eff_xtrl_tight                             = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_eff_xtrl_tight"));
        auto data_nbarlayer13_eff_xtrl_tight                        = static_cast<TEfficiency*>(datafile->Get("efficiencies/nbarlayer13_eff_xtrl_tight"));
        auto data_maxrms_and_nbarlayer13_eff_xtrl_tight             = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_and_nbarlayer13_eff_xtrl_tight"));
        auto data_sumrms_low_energy_eff_xtrl_tight                  = static_cast<TEfficiency*>(datafile->Get("efficiencies/sumrms_low_energy_eff_xtrl_tight"));
        auto data_track_selection_eff_xtrl_tight                    = static_cast<TEfficiency*>(datafile->Get("efficiencies/track_selection_eff_xtrl_tight"));
        auto data_track_selection_eff_within_stk_fvolume_xtrl_tight = static_cast<TEfficiency*>(datafile->Get("efficiencies/track_selection_eff_within_stk_fvolume_xtrl_tight"));
        auto data_stk_1_rm_eff_xtrl_tight                           = static_cast<TEfficiency*>(datafile->Get("efficiencies/stk_1_rm_eff_xtrl_tight"));
        auto data_psd_stk_match_eff_xtrl_tight                      = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_stk_match_eff_xtrl_tight"));
        auto data_psd_charge_eff_xtrl_tight                         = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_charge_eff_xtrl_tight"));
        auto data_stk_charge_eff_xtrl_tight                         = static_cast<TEfficiency*>(datafile->Get("efficiencies/stk_charge_eff_xtrl_tight"));

        auto data_trigger_eff_het_over_let_xtrl_loose               = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_let_xtrl_loose"));
        auto data_trigger_eff_het_over_unb_xtrl_loose               = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_unb_xtrl_loose"));
        auto data_maxrms_eff_xtrl_loose                             = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_eff_xtrl_loose"));
        auto data_nbarlayer13_eff_xtrl_loose                        = static_cast<TEfficiency*>(datafile->Get("efficiencies/nbarlayer13_eff_xtrl_loose"));
        auto data_maxrms_and_nbarlayer13_eff_xtrl_loose             = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_and_nbarlayer13_eff_xtrl_loose"));
        auto data_sumrms_low_energy_eff_xtrl_loose                  = static_cast<TEfficiency*>(datafile->Get("efficiencies/sumrms_low_energy_eff_xtrl_loose"));
        auto data_track_selection_eff_xtrl_loose                    = static_cast<TEfficiency*>(datafile->Get("efficiencies/track_selection_eff_xtrl_loose"));
        auto data_track_selection_eff_within_stk_fvolume_xtrl_loose = static_cast<TEfficiency*>(datafile->Get("efficiencies/track_selection_eff_within_stk_fvolume_xtrl_loose"));
        auto data_stk_1_rm_eff_xtrl_loose                           = static_cast<TEfficiency*>(datafile->Get("efficiencies/stk_1_rm_eff_xtrl_loose"));
        auto data_psd_stk_match_eff_xtrl_loose                      = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_stk_match_eff_xtrl_loose"));
        auto data_psd_charge_eff_xtrl_loose                         = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_charge_eff_xtrl_loose"));
        auto data_stk_charge_eff_xtrl_loose                         = static_cast<TEfficiency*>(datafile->Get("efficiencies/stk_charge_eff_xtrl_loose"));

        auto data_trigger_eff_het_over_let_bdt                      = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_let_bdt"));
        auto data_trigger_eff_het_over_unb_bdt                      = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_unb_bdt"));
        auto data_maxrms_eff_bdt                                    = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_eff_bdt"));
        auto data_nbarlayer13_eff_bdt                               = static_cast<TEfficiency*>(datafile->Get("efficiencies/nbarlayer13_eff_bdt"));
        auto data_maxrms_and_nbarlayer13_eff_bdt                    = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_and_nbarlayer13_eff_bdt"));
        auto data_sumrms_low_energy_eff_bdt                         = static_cast<TEfficiency*>(datafile->Get("efficiencies/sumrms_low_energy_eff_bdt"));
        auto data_track_selection_eff_bdt                           = static_cast<TEfficiency*>(datafile->Get("efficiencies/track_selection_eff_bdt"));
        auto data_track_selection_eff_within_stk_fvolume_bdt        = static_cast<TEfficiency*>(datafile->Get("efficiencies/track_selection_eff_within_stk_fvolume_bdt"));
        auto data_stk_1_rm_eff_bdt                                  = static_cast<TEfficiency*>(datafile->Get("efficiencies/stk_1_rm_eff_bdt"));
        auto data_psd_stk_match_eff_bdt                             = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_stk_match_eff_bdt"));
        auto data_psd_charge_eff_bdt                                = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_charge_eff_bdt"));
        auto data_stk_charge_eff_bdt                                = static_cast<TEfficiency*>(datafile->Get("efficiencies/stk_charge_eff_bdt"));

        auto data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight     = static_cast<TEfficiency*>(datafile->Get("efficiencies/clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight"));
        auto data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight    = static_cast<TEfficiency*>(datafile->Get("efficiencies/clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight"));
        auto data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose     = static_cast<TEfficiency*>(datafile->Get("efficiencies/clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose"));
        auto data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose    = static_cast<TEfficiency*>(datafile->Get("efficiencies/clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose"));
        auto data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt       = static_cast<TEfficiency*>(datafile->Get("efficiencies/clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt"));
        auto data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt      = static_cast<TEfficiency*>(datafile->Get("efficiencies/clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt"));

        data_trigger_eff_het_over_let_xtrl_tight                ->SetLineColor(kRed);
        data_trigger_eff_het_over_unb_xtrl_tight                ->SetLineColor(kRed);
        data_maxrms_eff_xtrl_tight                              ->SetLineColor(kRed);
        data_nbarlayer13_eff_xtrl_tight                         ->SetLineColor(kRed);
        data_maxrms_and_nbarlayer13_eff_xtrl_tight              ->SetLineColor(kRed);
        data_sumrms_low_energy_eff_xtrl_tight                   ->SetLineColor(kRed);
        data_track_selection_eff_xtrl_tight                     ->SetLineColor(kRed);
        data_track_selection_eff_within_stk_fvolume_xtrl_tight  ->SetLineColor(kRed);
        data_stk_1_rm_eff_xtrl_tight                            ->SetLineColor(kRed);
        data_psd_stk_match_eff_xtrl_tight                       ->SetLineColor(kRed);
        data_psd_charge_eff_xtrl_tight                          ->SetLineColor(kRed);
        data_stk_charge_eff_xtrl_tight                          ->SetLineColor(kRed);

        data_trigger_eff_het_over_let_xtrl_loose                ->SetLineColor(kRed);
        data_trigger_eff_het_over_unb_xtrl_loose                ->SetLineColor(kRed);
        data_maxrms_eff_xtrl_loose                              ->SetLineColor(kRed);
        data_nbarlayer13_eff_xtrl_loose                         ->SetLineColor(kRed);
        data_maxrms_and_nbarlayer13_eff_xtrl_loose              ->SetLineColor(kRed);
        data_sumrms_low_energy_eff_xtrl_loose                   ->SetLineColor(kRed);
        data_track_selection_eff_xtrl_loose                     ->SetLineColor(kRed);
        data_track_selection_eff_within_stk_fvolume_xtrl_loose  ->SetLineColor(kRed);
        data_stk_1_rm_eff_xtrl_loose                            ->SetLineColor(kRed);
        data_psd_stk_match_eff_xtrl_loose                       ->SetLineColor(kRed);
        data_psd_charge_eff_xtrl_loose                          ->SetLineColor(kRed);
        data_stk_charge_eff_xtrl_loose                          ->SetLineColor(kRed);

        data_trigger_eff_het_over_let_bdt                       ->SetLineColor(kMagenta);
        data_trigger_eff_het_over_unb_bdt                       ->SetLineColor(kMagenta);
        data_maxrms_eff_bdt                                     ->SetLineColor(kMagenta);
        data_nbarlayer13_eff_bdt                                ->SetLineColor(kMagenta);
        data_maxrms_and_nbarlayer13_eff_bdt                     ->SetLineColor(kMagenta);
        data_sumrms_low_energy_eff_bdt                          ->SetLineColor(kMagenta);
        data_track_selection_eff_bdt                            ->SetLineColor(kMagenta);
        data_track_selection_eff_within_stk_fvolume_bdt         ->SetLineColor(kMagenta);
        data_stk_1_rm_eff_bdt                                   ->SetLineColor(kMagenta);
        data_psd_stk_match_eff_bdt                              ->SetLineColor(kMagenta);
        data_psd_charge_eff_bdt                                 ->SetLineColor(kMagenta);
        data_stk_charge_eff_bdt                                 ->SetLineColor(kMagenta);

        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight      ->SetLineColor(kRed);
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight     ->SetLineColor(kRed);
        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose      ->SetLineColor(kRed);
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose     ->SetLineColor(kRed);
        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt        ->SetLineColor(kMagenta);
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt       ->SetLineColor(kMagenta);

        data_trigger_eff_het_over_let_xtrl_tight                ->SetLineWidth(2);
        data_trigger_eff_het_over_unb_xtrl_tight                ->SetLineWidth(2);
        data_maxrms_eff_xtrl_tight                              ->SetLineWidth(2);
        data_nbarlayer13_eff_xtrl_tight                         ->SetLineWidth(2);
        data_maxrms_and_nbarlayer13_eff_xtrl_tight              ->SetLineWidth(2);
        data_sumrms_low_energy_eff_xtrl_tight                   ->SetLineWidth(2);
        data_track_selection_eff_xtrl_tight                     ->SetLineWidth(2);
        data_track_selection_eff_within_stk_fvolume_xtrl_tight  ->SetLineWidth(2);
        data_stk_1_rm_eff_xtrl_tight                            ->SetLineWidth(2);
        data_psd_stk_match_eff_xtrl_tight                       ->SetLineWidth(2);
        data_psd_charge_eff_xtrl_tight                          ->SetLineWidth(2);
        data_stk_charge_eff_xtrl_tight                          ->SetLineWidth(2);

        data_trigger_eff_het_over_let_xtrl_loose                ->SetLineWidth(2);
        data_trigger_eff_het_over_unb_xtrl_loose                ->SetLineWidth(2);
        data_maxrms_eff_xtrl_loose                              ->SetLineWidth(2);
        data_nbarlayer13_eff_xtrl_loose                         ->SetLineWidth(2);
        data_maxrms_and_nbarlayer13_eff_xtrl_loose              ->SetLineWidth(2);
        data_sumrms_low_energy_eff_xtrl_loose                   ->SetLineWidth(2);
        data_track_selection_eff_xtrl_loose                     ->SetLineWidth(2);
        data_track_selection_eff_within_stk_fvolume_xtrl_loose  ->SetLineWidth(2);
        data_stk_1_rm_eff_xtrl_loose                            ->SetLineWidth(2);
        data_psd_stk_match_eff_xtrl_loose                       ->SetLineWidth(2);
        data_psd_charge_eff_xtrl_loose                          ->SetLineWidth(2);
        data_stk_charge_eff_xtrl_loose                          ->SetLineWidth(2);

        data_trigger_eff_het_over_let_bdt                       ->SetLineWidth(2);
        data_trigger_eff_het_over_unb_bdt                       ->SetLineWidth(2);
        data_maxrms_eff_bdt                                     ->SetLineWidth(2);
        data_nbarlayer13_eff_bdt                                ->SetLineWidth(2);
        data_maxrms_and_nbarlayer13_eff_bdt                     ->SetLineWidth(2);
        data_sumrms_low_energy_eff_bdt                          ->SetLineWidth(2);
        data_track_selection_eff_bdt                            ->SetLineWidth(2);
        data_track_selection_eff_within_stk_fvolume_bdt         ->SetLineWidth(2);
        data_stk_1_rm_eff_bdt                                   ->SetLineWidth(2);
        data_psd_stk_match_eff_bdt                              ->SetLineWidth(2);
        data_psd_charge_eff_bdt                                 ->SetLineWidth(2);
        data_stk_charge_eff_bdt                                 ->SetLineWidth(2);

        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight      ->SetLineWidth(2);
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight     ->SetLineWidth(2);
        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose      ->SetLineWidth(2);
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose     ->SetLineWidth(2);
        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt        ->SetLineWidth(2);
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt       ->SetLineWidth(2);

        data_trigger_eff_het_over_let_xtrl_tight                ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_trigger_eff_het_over_unb_xtrl_tight                ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_maxrms_eff_xtrl_tight                              ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_nbarlayer13_eff_xtrl_tight                         ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_maxrms_and_nbarlayer13_eff_xtrl_tight              ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_sumrms_low_energy_eff_xtrl_tight                   ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_track_selection_eff_xtrl_tight                     ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_track_selection_eff_within_stk_fvolume_xtrl_tight  ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_stk_1_rm_eff_xtrl_tight                            ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_psd_stk_match_eff_xtrl_tight                       ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_psd_charge_eff_xtrl_tight                          ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_stk_charge_eff_xtrl_tight                          ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());

        data_trigger_eff_het_over_let_xtrl_loose                ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_trigger_eff_het_over_unb_xtrl_loose                ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_maxrms_eff_xtrl_loose                              ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_nbarlayer13_eff_xtrl_loose                         ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_maxrms_and_nbarlayer13_eff_xtrl_loose              ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_sumrms_low_energy_eff_xtrl_loose                   ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_track_selection_eff_xtrl_loose                     ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_track_selection_eff_within_stk_fvolume_xtrl_loose  ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_stk_1_rm_eff_xtrl_loose                            ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_psd_stk_match_eff_xtrl_loose                       ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_psd_charge_eff_xtrl_loose                          ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_stk_charge_eff_xtrl_loose                          ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());

        data_trigger_eff_het_over_let_bdt                       ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_trigger_eff_het_over_unb_bdt                       ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_maxrms_eff_bdt                                     ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_nbarlayer13_eff_bdt                                ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_maxrms_and_nbarlayer13_eff_bdt                     ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_sumrms_low_energy_eff_bdt                          ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_track_selection_eff_bdt                            ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_track_selection_eff_within_stk_fvolume_bdt         ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_stk_1_rm_eff_bdt                                   ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_psd_stk_match_eff_bdt                              ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_psd_charge_eff_bdt                                 ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_stk_charge_eff_bdt                                 ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());

        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight      ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight     ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose      ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose     ->SetTitle((std::string(data_label) + std::string(" (xtrl)")).c_str());
        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt        ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt       ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());

        data_trigger_eff_het_over_let_xtrl_tight                ->SetDirectory(0);
        data_trigger_eff_het_over_unb_xtrl_tight                ->SetDirectory(0);
        data_maxrms_eff_xtrl_tight                              ->SetDirectory(0);
        data_nbarlayer13_eff_xtrl_tight                         ->SetDirectory(0);
        data_maxrms_and_nbarlayer13_eff_xtrl_tight              ->SetDirectory(0);
        data_sumrms_low_energy_eff_xtrl_tight                   ->SetDirectory(0);
        data_track_selection_eff_xtrl_tight                     ->SetDirectory(0);
        data_track_selection_eff_within_stk_fvolume_xtrl_tight  ->SetDirectory(0);
        data_stk_1_rm_eff_xtrl_tight                            ->SetDirectory(0);
        data_psd_stk_match_eff_xtrl_tight                       ->SetDirectory(0);
        data_psd_charge_eff_xtrl_tight                          ->SetDirectory(0);
        data_stk_charge_eff_xtrl_tight                          ->SetDirectory(0);

        data_trigger_eff_het_over_let_xtrl_loose                ->SetDirectory(0);
        data_trigger_eff_het_over_unb_xtrl_loose                ->SetDirectory(0);
        data_maxrms_eff_xtrl_loose                              ->SetDirectory(0);
        data_nbarlayer13_eff_xtrl_loose                         ->SetDirectory(0);
        data_maxrms_and_nbarlayer13_eff_xtrl_loose              ->SetDirectory(0);
        data_sumrms_low_energy_eff_xtrl_loose                   ->SetDirectory(0);
        data_track_selection_eff_xtrl_loose                     ->SetDirectory(0);
        data_track_selection_eff_within_stk_fvolume_xtrl_loose  ->SetDirectory(0);
        data_stk_1_rm_eff_xtrl_loose                            ->SetDirectory(0);
        data_psd_stk_match_eff_xtrl_loose                       ->SetDirectory(0);
        data_psd_charge_eff_xtrl_loose                          ->SetDirectory(0);
        data_stk_charge_eff_xtrl_loose                          ->SetDirectory(0);

        data_trigger_eff_het_over_let_bdt                       ->SetDirectory(0);
        data_trigger_eff_het_over_unb_bdt                       ->SetDirectory(0);
        data_maxrms_eff_bdt                                     ->SetDirectory(0);
        data_nbarlayer13_eff_bdt                                ->SetDirectory(0);
        data_maxrms_and_nbarlayer13_eff_bdt                     ->SetDirectory(0);
        data_sumrms_low_energy_eff_bdt                          ->SetDirectory(0);
        data_track_selection_eff_bdt                            ->SetDirectory(0);
        data_track_selection_eff_within_stk_fvolume_bdt         ->SetDirectory(0);
        data_stk_1_rm_eff_bdt                                   ->SetDirectory(0);
        data_psd_stk_match_eff_bdt                              ->SetDirectory(0);
        data_psd_charge_eff_bdt                                 ->SetDirectory(0);
        data_stk_charge_eff_bdt                                 ->SetDirectory(0);

        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight      ->SetDirectory(0);
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight     ->SetDirectory(0);
        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose      ->SetDirectory(0);
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose     ->SetDirectory(0);
        data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt        ->SetDirectory(0);
        data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt       ->SetDirectory(0);

        datafile->Close();

        TFile* mcfile = TFile::Open(mc_file, "READ");
        if (mcfile->IsZombie()) {
            std::cerr << "\n\nError opening input MC ROOT file [" << mc_file << "]\n\n";
            exit(100);
        }

        auto mc_trigger_eff_het_over_let_xtrl_tight                 = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_let_xtrl_tight"));
        auto mc_trigger_eff_het_over_unb_xtrl_tight                 = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_unb_xtrl_tight"));
        auto mc_maxrms_eff_xtrl_tight                               = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_eff_xtrl_tight"));
        auto mc_nbarlayer13_eff_xtrl_tight                          = static_cast<TEfficiency*>(mcfile->Get("efficiencies/nbarlayer13_eff_xtrl_tight"));
        auto mc_maxrms_and_nbarlayer13_eff_xtrl_tight               = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_and_nbarlayer13_eff_xtrl_tight"));
        auto mc_sumrms_low_energy_eff_xtrl_tight                    = static_cast<TEfficiency*>(mcfile->Get("efficiencies/sumrms_low_energy_eff_xtrl_tight"));
        auto mc_track_selection_eff_xtrl_tight                      = static_cast<TEfficiency*>(mcfile->Get("efficiencies/track_selection_eff_xtrl_tight"));
        auto mc_track_selection_eff_within_stk_fvolume_xtrl_tight   = static_cast<TEfficiency*>(mcfile->Get("efficiencies/track_selection_eff_within_stk_fvolume_xtrl_tight"));
        auto mc_stk_1_rm_eff_xtrl_tight                             = static_cast<TEfficiency*>(mcfile->Get("efficiencies/stk_1_rm_eff_xtrl_tight"));
        auto mc_psd_stk_match_eff_xtrl_tight                        = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_stk_match_eff_xtrl_tight"));
        auto mc_psd_charge_eff_xtrl_tight                           = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_charge_eff_xtrl_tight"));
        auto mc_stk_charge_eff_xtrl_tight                           = static_cast<TEfficiency*>(mcfile->Get("efficiencies/stk_charge_eff_xtrl_tight"));

        auto mc_trigger_eff_het_over_let_xtrl_loose                 = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_let_xtrl_loose"));
        auto mc_trigger_eff_het_over_unb_xtrl_loose                 = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_unb_xtrl_loose"));
        auto mc_maxrms_eff_xtrl_loose                               = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_eff_xtrl_loose"));
        auto mc_nbarlayer13_eff_xtrl_loose                          = static_cast<TEfficiency*>(mcfile->Get("efficiencies/nbarlayer13_eff_xtrl_loose"));
        auto mc_maxrms_and_nbarlayer13_eff_xtrl_loose               = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_and_nbarlayer13_eff_xtrl_loose"));
        auto mc_sumrms_low_energy_eff_xtrl_loose                    = static_cast<TEfficiency*>(mcfile->Get("efficiencies/sumrms_low_energy_eff_xtrl_loose"));
        auto mc_track_selection_eff_xtrl_loose                      = static_cast<TEfficiency*>(mcfile->Get("efficiencies/track_selection_eff_xtrl_loose"));
        auto mc_track_selection_eff_within_stk_fvolume_xtrl_loose   = static_cast<TEfficiency*>(mcfile->Get("efficiencies/track_selection_eff_within_stk_fvolume_xtrl_loose"));
        auto mc_stk_1_rm_eff_xtrl_loose                             = static_cast<TEfficiency*>(mcfile->Get("efficiencies/stk_1_rm_eff_xtrl_loose"));
        auto mc_psd_stk_match_eff_xtrl_loose                        = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_stk_match_eff_xtrl_loose"));
        auto mc_psd_charge_eff_xtrl_loose                           = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_charge_eff_xtrl_loose"));
        auto mc_stk_charge_eff_xtrl_loose                           = static_cast<TEfficiency*>(mcfile->Get("efficiencies/stk_charge_eff_xtrl_loose"));
        
        auto mc_trigger_eff_het_over_let_bdt                        = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_let_bdt"));
        auto mc_trigger_eff_het_over_unb_bdt                        = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_unb_bdt"));
        auto mc_maxrms_eff_bdt                                      = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_eff_bdt"));
        auto mc_nbarlayer13_eff_bdt                                 = static_cast<TEfficiency*>(mcfile->Get("efficiencies/nbarlayer13_eff_bdt"));
        auto mc_maxrms_and_nbarlayer13_eff_bdt                      = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_and_nbarlayer13_eff_bdt"));
        auto mc_sumrms_low_energy_eff_bdt                           = static_cast<TEfficiency*>(mcfile->Get("efficiencies/sumrms_low_energy_eff_bdt"));
        auto mc_track_selection_eff_bdt                             = static_cast<TEfficiency*>(mcfile->Get("efficiencies/track_selection_eff_bdt"));
        auto mc_track_selection_eff_within_stk_fvolume_bdt          = static_cast<TEfficiency*>(mcfile->Get("efficiencies/track_selection_eff_within_stk_fvolume_bdt"));
        auto mc_stk_1_rm_eff_bdt                                    = static_cast<TEfficiency*>(mcfile->Get("efficiencies/stk_1_rm_eff_bdt"));
        auto mc_psd_stk_match_eff_bdt                               = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_stk_match_eff_bdt"));
        auto mc_psd_charge_eff_bdt                                  = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_charge_eff_bdt"));
        auto mc_stk_charge_eff_bdt                                  = static_cast<TEfficiency*>(mcfile->Get("efficiencies/stk_charge_eff_bdt"));

        auto mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight     = static_cast<TEfficiency*>(mcfile->Get("efficiencies/clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight"));
        auto mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight    = static_cast<TEfficiency*>(mcfile->Get("efficiencies/clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight"));
        auto mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose     = static_cast<TEfficiency*>(mcfile->Get("efficiencies/clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose"));
        auto mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose    = static_cast<TEfficiency*>(mcfile->Get("efficiencies/clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose"));
        auto mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt       = static_cast<TEfficiency*>(mcfile->Get("efficiencies/clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt"));
        auto mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt      = static_cast<TEfficiency*>(mcfile->Get("efficiencies/clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt"));

        mc_trigger_eff_het_over_let_xtrl_tight                  ->SetLineColor(kBlue);
        mc_trigger_eff_het_over_unb_xtrl_tight                  ->SetLineColor(kBlue);
        mc_maxrms_eff_xtrl_tight                                ->SetLineColor(kBlue);
        mc_nbarlayer13_eff_xtrl_tight                           ->SetLineColor(kBlue);
        mc_maxrms_and_nbarlayer13_eff_xtrl_tight                ->SetLineColor(kBlue);
        mc_sumrms_low_energy_eff_xtrl_tight                     ->SetLineColor(kBlue);
        mc_track_selection_eff_xtrl_tight                       ->SetLineColor(kBlue);
        mc_track_selection_eff_within_stk_fvolume_xtrl_tight    ->SetLineColor(kBlue);
        mc_stk_1_rm_eff_xtrl_tight                              ->SetLineColor(kBlue);
        mc_psd_stk_match_eff_xtrl_tight                         ->SetLineColor(kBlue);
        mc_psd_charge_eff_xtrl_tight                            ->SetLineColor(kBlue);
        mc_stk_charge_eff_xtrl_tight                            ->SetLineColor(kBlue);

        mc_trigger_eff_het_over_let_xtrl_loose                  ->SetLineColor(kBlue);
        mc_trigger_eff_het_over_unb_xtrl_loose                  ->SetLineColor(kBlue);
        mc_maxrms_eff_xtrl_loose                                ->SetLineColor(kBlue);
        mc_nbarlayer13_eff_xtrl_loose                           ->SetLineColor(kBlue);
        mc_maxrms_and_nbarlayer13_eff_xtrl_loose                ->SetLineColor(kBlue);
        mc_sumrms_low_energy_eff_xtrl_loose                     ->SetLineColor(kBlue);
        mc_track_selection_eff_xtrl_loose                       ->SetLineColor(kBlue);
        mc_track_selection_eff_within_stk_fvolume_xtrl_loose    ->SetLineColor(kBlue);
        mc_stk_1_rm_eff_xtrl_loose                              ->SetLineColor(kBlue);
        mc_psd_stk_match_eff_xtrl_loose                         ->SetLineColor(kBlue);
        mc_psd_charge_eff_xtrl_loose                            ->SetLineColor(kBlue);
        mc_stk_charge_eff_xtrl_loose                            ->SetLineColor(kBlue);

        mc_trigger_eff_het_over_let_bdt                         ->SetLineColor(kGreen);
        mc_trigger_eff_het_over_unb_bdt                         ->SetLineColor(kGreen);
        mc_maxrms_eff_bdt                                       ->SetLineColor(kGreen);
        mc_nbarlayer13_eff_bdt                                  ->SetLineColor(kGreen);
        mc_maxrms_and_nbarlayer13_eff_bdt                       ->SetLineColor(kGreen);
        mc_sumrms_low_energy_eff_bdt                            ->SetLineColor(kGreen);
        mc_track_selection_eff_bdt                              ->SetLineColor(kGreen);
        mc_track_selection_eff_within_stk_fvolume_bdt           ->SetLineColor(kGreen);
        mc_stk_1_rm_eff_bdt                                     ->SetLineColor(kGreen);
        mc_psd_stk_match_eff_bdt                                ->SetLineColor(kGreen);
        mc_psd_charge_eff_bdt                                   ->SetLineColor(kGreen);
        mc_stk_charge_eff_bdt                                   ->SetLineColor(kGreen);

        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight      ->SetLineColor(kBlue);
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight     ->SetLineColor(kBlue);
        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose      ->SetLineColor(kBlue);
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose     ->SetLineColor(kBlue);
        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt        ->SetLineColor(kGreen);
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt       ->SetLineColor(kGreen);

        mc_trigger_eff_het_over_let_xtrl_tight                  ->SetLineWidth(2);
        mc_trigger_eff_het_over_unb_xtrl_tight                  ->SetLineWidth(2);
        mc_maxrms_eff_xtrl_tight                                ->SetLineWidth(2);
        mc_nbarlayer13_eff_xtrl_tight                           ->SetLineWidth(2);
        mc_maxrms_and_nbarlayer13_eff_xtrl_tight                ->SetLineWidth(2);
        mc_sumrms_low_energy_eff_xtrl_tight                     ->SetLineWidth(2);
        mc_track_selection_eff_xtrl_tight                       ->SetLineWidth(2);
        mc_track_selection_eff_within_stk_fvolume_xtrl_tight    ->SetLineWidth(2);
        mc_stk_1_rm_eff_xtrl_tight                              ->SetLineWidth(2);
        mc_psd_stk_match_eff_xtrl_tight                         ->SetLineWidth(2);
        mc_psd_charge_eff_xtrl_tight                            ->SetLineWidth(2);
        mc_stk_charge_eff_xtrl_tight                            ->SetLineWidth(2);

        mc_trigger_eff_het_over_let_xtrl_loose                  ->SetLineWidth(2);
        mc_trigger_eff_het_over_unb_xtrl_loose                  ->SetLineWidth(2);
        mc_maxrms_eff_xtrl_loose                                ->SetLineWidth(2);
        mc_nbarlayer13_eff_xtrl_loose                           ->SetLineWidth(2);
        mc_maxrms_and_nbarlayer13_eff_xtrl_loose                ->SetLineWidth(2);
        mc_sumrms_low_energy_eff_xtrl_loose                     ->SetLineWidth(2);
        mc_track_selection_eff_xtrl_loose                       ->SetLineWidth(2);
        mc_track_selection_eff_within_stk_fvolume_xtrl_loose    ->SetLineWidth(2);
        mc_stk_1_rm_eff_xtrl_loose                              ->SetLineWidth(2);
        mc_psd_stk_match_eff_xtrl_loose                         ->SetLineWidth(2);
        mc_psd_charge_eff_xtrl_loose                            ->SetLineWidth(2);
        mc_stk_charge_eff_xtrl_loose                            ->SetLineWidth(2);

        mc_trigger_eff_het_over_let_bdt                         ->SetLineWidth(2);
        mc_trigger_eff_het_over_unb_bdt                         ->SetLineWidth(2);
        mc_maxrms_eff_bdt                                       ->SetLineWidth(2);
        mc_nbarlayer13_eff_bdt                                  ->SetLineWidth(2);
        mc_maxrms_and_nbarlayer13_eff_bdt                       ->SetLineWidth(2);
        mc_sumrms_low_energy_eff_bdt                            ->SetLineWidth(2);
        mc_track_selection_eff_bdt                              ->SetLineWidth(2);
        mc_track_selection_eff_within_stk_fvolume_bdt           ->SetLineWidth(2);
        mc_stk_1_rm_eff_bdt                                     ->SetLineWidth(2);
        mc_psd_stk_match_eff_bdt                                ->SetLineWidth(2);
        mc_psd_charge_eff_bdt                                   ->SetLineWidth(2);
        mc_stk_charge_eff_bdt                                   ->SetLineWidth(2);

        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight      ->SetLineWidth(2);
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight     ->SetLineWidth(2);
        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose      ->SetLineWidth(2);
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose     ->SetLineWidth(2);
        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt        ->SetLineWidth(2);
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt       ->SetLineWidth(2);

        mc_trigger_eff_het_over_let_xtrl_tight                  ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_trigger_eff_het_over_unb_xtrl_tight                  ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_maxrms_eff_xtrl_tight                                ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_nbarlayer13_eff_xtrl_tight                           ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_maxrms_and_nbarlayer13_eff_xtrl_tight                ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_sumrms_low_energy_eff_xtrl_tight                     ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_track_selection_eff_xtrl_tight                       ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_track_selection_eff_within_stk_fvolume_xtrl_tight    ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_stk_1_rm_eff_xtrl_tight                              ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_psd_stk_match_eff_xtrl_tight                         ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_psd_charge_eff_xtrl_tight                            ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_stk_charge_eff_xtrl_tight                            ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());

        mc_trigger_eff_het_over_let_xtrl_loose                  ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_trigger_eff_het_over_unb_xtrl_loose                  ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_maxrms_eff_xtrl_loose                                ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_nbarlayer13_eff_xtrl_loose                           ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_maxrms_and_nbarlayer13_eff_xtrl_loose                ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_sumrms_low_energy_eff_xtrl_loose                     ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_track_selection_eff_xtrl_loose                       ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_track_selection_eff_within_stk_fvolume_xtrl_loose    ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_stk_1_rm_eff_xtrl_loose                              ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_psd_stk_match_eff_xtrl_loose                         ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_psd_charge_eff_xtrl_loose                            ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_stk_charge_eff_xtrl_loose                            ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());

        mc_trigger_eff_het_over_let_bdt                         ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_trigger_eff_het_over_unb_bdt                         ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_maxrms_eff_bdt                                       ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_nbarlayer13_eff_bdt                                  ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_maxrms_and_nbarlayer13_eff_bdt                       ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_sumrms_low_energy_eff_bdt                            ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_track_selection_eff_bdt                              ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_track_selection_eff_within_stk_fvolume_bdt           ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_stk_1_rm_eff_bdt                                     ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_psd_stk_match_eff_bdt                                ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_psd_charge_eff_bdt                                   ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_stk_charge_eff_bdt                                   ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());

        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight      ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight     ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose      ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose     ->SetTitle((std::string(mc_label) + std::string(" (xtrl)")).c_str());
        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt        ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt       ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());

        mc_trigger_eff_het_over_let_xtrl_tight                  ->SetDirectory(0);
        mc_trigger_eff_het_over_unb_xtrl_tight                  ->SetDirectory(0);
        mc_maxrms_eff_xtrl_tight                                ->SetDirectory(0);
        mc_nbarlayer13_eff_xtrl_tight                           ->SetDirectory(0);
        mc_maxrms_and_nbarlayer13_eff_xtrl_tight                ->SetDirectory(0);
        mc_sumrms_low_energy_eff_xtrl_tight                     ->SetDirectory(0);
        mc_track_selection_eff_xtrl_tight                       ->SetDirectory(0);
        mc_track_selection_eff_within_stk_fvolume_xtrl_tight    ->SetDirectory(0);
        mc_stk_1_rm_eff_xtrl_tight                              ->SetDirectory(0);
        mc_psd_stk_match_eff_xtrl_tight                         ->SetDirectory(0);
        mc_psd_charge_eff_xtrl_tight                            ->SetDirectory(0);
        mc_stk_charge_eff_xtrl_tight                            ->SetDirectory(0);

        mc_trigger_eff_het_over_let_xtrl_loose                  ->SetDirectory(0);
        mc_trigger_eff_het_over_unb_xtrl_loose                  ->SetDirectory(0);
        mc_maxrms_eff_xtrl_loose                                ->SetDirectory(0);
        mc_nbarlayer13_eff_xtrl_loose                           ->SetDirectory(0);
        mc_maxrms_and_nbarlayer13_eff_xtrl_loose                ->SetDirectory(0);
        mc_sumrms_low_energy_eff_xtrl_loose                     ->SetDirectory(0);
        mc_track_selection_eff_xtrl_loose                       ->SetDirectory(0);
        mc_track_selection_eff_within_stk_fvolume_xtrl_loose    ->SetDirectory(0);
        mc_stk_1_rm_eff_xtrl_loose                              ->SetDirectory(0);
        mc_psd_stk_match_eff_xtrl_loose                         ->SetDirectory(0);
        mc_psd_charge_eff_xtrl_loose                            ->SetDirectory(0);
        mc_stk_charge_eff_xtrl_loose                            ->SetDirectory(0);

        mc_trigger_eff_het_over_let_bdt                         ->SetDirectory(0);
        mc_trigger_eff_het_over_unb_bdt                         ->SetDirectory(0);
        mc_maxrms_eff_bdt                                       ->SetDirectory(0);
        mc_nbarlayer13_eff_bdt                                  ->SetDirectory(0);
        mc_maxrms_and_nbarlayer13_eff_bdt                       ->SetDirectory(0);
        mc_sumrms_low_energy_eff_bdt                            ->SetDirectory(0);
        mc_track_selection_eff_bdt                              ->SetDirectory(0);
        mc_track_selection_eff_within_stk_fvolume_bdt           ->SetDirectory(0);
        mc_stk_1_rm_eff_bdt                                     ->SetDirectory(0);
        mc_psd_stk_match_eff_bdt                                ->SetDirectory(0);
        mc_psd_charge_eff_bdt                                   ->SetDirectory(0);
        mc_stk_charge_eff_bdt                                   ->SetDirectory(0);

        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight      ->SetDirectory(0);
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight     ->SetDirectory(0);
        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose      ->SetDirectory(0);
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose     ->SetDirectory(0);
        mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt        ->SetDirectory(0);
        mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt       ->SetDirectory(0);

        mcfile->Close();

        auto remove_x_error_from_TEfficiency = [](TEfficiency* eff) -> std::shared_ptr<TGraphAsymmErrors> {
            auto gr = std::make_shared<TGraphAsymmErrors>(*static_cast<TGraphAsymmErrors*>(eff->CreateGraph()));
            for (int i = 0; i < gr->GetN(); ++i) {
                gr->GetEXlow()[i] = 0;  
                gr->GetEXhigh()[i] = 0;  
            }
            return gr;
        };
        
        TCanvas print_canvas("print_canvas", "print_canvas");
        print_canvas.SetTicks();

        // Print trigger HET over LET efficiency (xtrl tight)
        auto gr_data_trigger_eff_het_over_let_xtrl_tight = remove_x_error_from_TEfficiency(data_trigger_eff_het_over_let_xtrl_tight);
        auto gr_mc_trigger_eff_het_over_let_xtrl_tight = remove_x_error_from_TEfficiency(mc_trigger_eff_het_over_let_xtrl_tight);
        auto gr_data_trigger_eff_het_over_let_bdt = remove_x_error_from_TEfficiency(data_trigger_eff_het_over_let_bdt);
        auto gr_mc_trigger_eff_het_over_let_bdt = remove_x_error_from_TEfficiency(mc_trigger_eff_het_over_let_bdt);

        gr_data_trigger_eff_het_over_let_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_trigger_eff_het_over_let_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_trigger_eff_het_over_let_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_trigger_eff_het_over_let_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_trigger_eff_het_over_let_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_trigger_eff_het_over_let_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_trigger_eff_het_over_let_bdt->SetMarkerColor(kMagenta);
        gr_mc_trigger_eff_het_over_let_bdt->SetMarkerColor(kGreen);

        gr_data_trigger_eff_het_over_let_xtrl_tight->Draw("AP");
        gr_mc_trigger_eff_het_over_let_xtrl_tight->Draw("P");
        gr_data_trigger_eff_het_over_let_bdt->Draw("P");
        gr_mc_trigger_eff_het_over_let_bdt->Draw("P");

        gPad->Update(); 
        gr_data_trigger_eff_het_over_let_xtrl_tight->SetMinimum(0.5);
        gr_data_trigger_eff_het_over_let_xtrl_tight->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        auto legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        TPaveLabel label(0.0, 0.95, 0.3, 1, "Trigger Efficiency - HET over LET (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf(","Title:Trigger HET over LET - xtrl tight efficiency");
        
        // Print trigger HET over LET efficiency (xtrl loose)
        auto gr_data_trigger_eff_het_over_let_xtrl_loose = remove_x_error_from_TEfficiency(data_trigger_eff_het_over_let_xtrl_loose);
        auto gr_mc_trigger_eff_het_over_let_xtrl_loose = remove_x_error_from_TEfficiency(mc_trigger_eff_het_over_let_xtrl_loose);

        gr_data_trigger_eff_het_over_let_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_trigger_eff_het_over_let_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_trigger_eff_het_over_let_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_trigger_eff_het_over_let_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_trigger_eff_het_over_let_xtrl_loose->Draw("AP");
        gr_mc_trigger_eff_het_over_let_xtrl_loose->Draw("P");
        gr_data_trigger_eff_het_over_let_bdt->Draw("P");
        gr_mc_trigger_eff_het_over_let_bdt->Draw("P");

        gPad->Update(); 
        gr_data_trigger_eff_het_over_let_xtrl_loose->SetMinimum(0.5);
        gr_data_trigger_eff_het_over_let_xtrl_loose->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "Trigger Efficiency - HET over LET (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:Trigger HET over LET - xtrl loose efficiency");
        
        // Print trigger HET over UNB efficiency (xtrl tight)
        auto gr_data_trigger_eff_het_over_unb_xtrl_tight = remove_x_error_from_TEfficiency(data_trigger_eff_het_over_unb_xtrl_tight);
        auto gr_mc_trigger_eff_het_over_unb_xtrl_tight = remove_x_error_from_TEfficiency(mc_trigger_eff_het_over_unb_xtrl_tight);
        auto gr_data_trigger_eff_het_over_unb_bdt = remove_x_error_from_TEfficiency(data_trigger_eff_het_over_unb_bdt);
        auto gr_mc_trigger_eff_het_over_unb_bdt = remove_x_error_from_TEfficiency(mc_trigger_eff_het_over_unb_bdt);

        gr_data_trigger_eff_het_over_unb_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_trigger_eff_het_over_unb_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_trigger_eff_het_over_unb_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_trigger_eff_het_over_unb_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_trigger_eff_het_over_unb_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_trigger_eff_het_over_unb_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_trigger_eff_het_over_unb_bdt->SetMarkerColor(kMagenta);
        gr_mc_trigger_eff_het_over_unb_bdt->SetMarkerColor(kGreen);

        gr_data_trigger_eff_het_over_unb_xtrl_tight->Draw("AP");
        gr_mc_trigger_eff_het_over_unb_xtrl_tight->Draw("P");
        gr_data_trigger_eff_het_over_unb_bdt->Draw("P");
        gr_mc_trigger_eff_het_over_unb_bdt->Draw("P");

        gPad->Update(); 
        gr_data_trigger_eff_het_over_unb_xtrl_tight->SetMinimum(0.5);
        gr_data_trigger_eff_het_over_unb_xtrl_tight->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "Trigger Efficiency - HET over UNB (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf","Title:Trigger HET over UNB - xtrl tight efficiency");
        
        // Print trigger HET over UNB efficiency (xtrl loose)
        auto gr_data_trigger_eff_het_over_unb_xtrl_loose = remove_x_error_from_TEfficiency(data_trigger_eff_het_over_unb_xtrl_loose);
        auto gr_mc_trigger_eff_het_over_unb_xtrl_loose = remove_x_error_from_TEfficiency(mc_trigger_eff_het_over_unb_xtrl_loose);

        gr_data_trigger_eff_het_over_unb_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_trigger_eff_het_over_unb_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_trigger_eff_het_over_unb_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_trigger_eff_het_over_unb_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_trigger_eff_het_over_unb_xtrl_loose->Draw("AP");
        gr_mc_trigger_eff_het_over_unb_xtrl_loose->Draw("P");
        gr_data_trigger_eff_het_over_unb_bdt->Draw("P");
        gr_mc_trigger_eff_het_over_unb_bdt->Draw("P");

        gPad->Update(); 
        gr_data_trigger_eff_het_over_unb_xtrl_tight->SetMinimum(0.5);
        gr_data_trigger_eff_het_over_unb_xtrl_tight->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "Trigger Efficiency - HET over UNB (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:Trigger HET over UNB - xtrl loose efficiency");

        // Print max RMS efficiency (xtrl tight)
        auto gr_data_maxrms_eff_xtrl_tight = remove_x_error_from_TEfficiency(data_maxrms_eff_xtrl_tight);
        auto gr_mc_maxrms_eff_xtrl_tight = remove_x_error_from_TEfficiency(mc_maxrms_eff_xtrl_tight);
        auto gr_data_maxrms_eff_bdt = remove_x_error_from_TEfficiency(data_maxrms_eff_bdt);
        auto gr_mc_maxrms_eff_bdt = remove_x_error_from_TEfficiency(mc_maxrms_eff_bdt);
        
        gr_data_maxrms_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_maxrms_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_maxrms_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_maxrms_eff_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_maxrms_eff_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_maxrms_eff_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_maxrms_eff_bdt->SetMarkerColor(kMagenta);
        gr_mc_maxrms_eff_bdt->SetMarkerColor(kGreen);

        gr_data_maxrms_eff_xtrl_tight->Draw("AP");
        gr_mc_maxrms_eff_xtrl_tight->Draw("P");
        gr_data_maxrms_eff_bdt->Draw("P");
        gr_mc_maxrms_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_maxrms_eff_xtrl_tight->SetMinimum(0.9);
        gr_data_maxrms_eff_xtrl_tight->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "maxRMS Efficiency (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf","Title:max RMS - xtrl tight efficiency");
        
        // Print max RMS efficiency (xtrl loose)
        auto gr_data_maxrms_eff_xtrl_loose = remove_x_error_from_TEfficiency(data_maxrms_eff_xtrl_loose);
        auto gr_mc_maxrms_eff_xtrl_loose = remove_x_error_from_TEfficiency(mc_maxrms_eff_xtrl_loose);

        gr_data_maxrms_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_maxrms_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_maxrms_eff_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_maxrms_eff_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_maxrms_eff_xtrl_loose->Draw("AP");
        gr_mc_maxrms_eff_xtrl_loose->Draw("P");
        gr_data_maxrms_eff_bdt->Draw("P");
        gr_mc_maxrms_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_maxrms_eff_xtrl_loose->SetMinimum(0.9);
        gr_data_maxrms_eff_xtrl_loose->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "maxRMS Efficiency (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:maxRMS - xtrl loose efficiency");

        // Print nbarlayer13 efficiency (xtrl tight)
        auto gr_data_nbarlayer13_eff_xtrl_tight = remove_x_error_from_TEfficiency(data_nbarlayer13_eff_xtrl_tight);
        auto gr_mc_nbarlayer13_eff_xtrl_tight = remove_x_error_from_TEfficiency(mc_nbarlayer13_eff_xtrl_tight);
        auto gr_data_nbarlayer13_eff_bdt = remove_x_error_from_TEfficiency(data_nbarlayer13_eff_bdt);
        auto gr_mc_nbarlayer13_eff_bdt = remove_x_error_from_TEfficiency(mc_nbarlayer13_eff_bdt);
        
        gr_data_nbarlayer13_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_nbarlayer13_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_nbarlayer13_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_nbarlayer13_eff_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_nbarlayer13_eff_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_nbarlayer13_eff_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_nbarlayer13_eff_bdt->SetMarkerColor(kMagenta);
        gr_mc_nbarlayer13_eff_bdt->SetMarkerColor(kGreen);

        gr_data_nbarlayer13_eff_xtrl_tight->Draw("AP");
        gr_mc_nbarlayer13_eff_xtrl_tight->Draw("P");
        gr_data_nbarlayer13_eff_bdt->Draw("P");
        gr_mc_nbarlayer13_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_nbarlayer13_eff_xtrl_tight->SetMinimum(0.9);
        gr_data_nbarlayer13_eff_xtrl_tight->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "nbarlayer13 Efficiency (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf","Title:nbarlayer13 - xtrl tight efficiency");
        
        // Print nbarlayer13 efficiency (xtrl loose)
        auto gr_data_nbarlayer13_eff_xtrl_loose = remove_x_error_from_TEfficiency(data_nbarlayer13_eff_xtrl_loose);
        auto gr_mc_nbarlayer13_eff_xtrl_loose = remove_x_error_from_TEfficiency(mc_nbarlayer13_eff_xtrl_loose);
        
        gr_data_nbarlayer13_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_nbarlayer13_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_nbarlayer13_eff_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_nbarlayer13_eff_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_nbarlayer13_eff_xtrl_loose->Draw("AP");
        gr_mc_nbarlayer13_eff_xtrl_loose->Draw("P");
        gr_data_nbarlayer13_eff_bdt->Draw("P");
        gr_mc_nbarlayer13_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_nbarlayer13_eff_xtrl_loose->SetMinimum(0.9);
        gr_data_nbarlayer13_eff_xtrl_loose->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "nbarlayer13 Efficiency (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:nbarlayer13 - xtrl loose efficiency");

        // Print maxrms and nbarlayer13 efficiency (xtrl tight)
        auto gr_data_maxrms_and_nbarlayer13_eff_xtrl_tight = remove_x_error_from_TEfficiency(data_maxrms_and_nbarlayer13_eff_xtrl_tight);
        auto gr_mc_maxrms_and_nbarlayer13_eff_xtrl_tight = remove_x_error_from_TEfficiency(mc_maxrms_and_nbarlayer13_eff_xtrl_tight);
        auto gr_data_maxrms_and_nbarlayer13_eff_bdt = remove_x_error_from_TEfficiency(data_maxrms_and_nbarlayer13_eff_bdt);
        auto gr_mc_maxrms_and_nbarlayer13_eff_bdt = remove_x_error_from_TEfficiency(mc_maxrms_and_nbarlayer13_eff_bdt);
        
        gr_data_maxrms_and_nbarlayer13_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_maxrms_and_nbarlayer13_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_maxrms_and_nbarlayer13_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_maxrms_and_nbarlayer13_eff_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_maxrms_and_nbarlayer13_eff_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_maxrms_and_nbarlayer13_eff_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_maxrms_and_nbarlayer13_eff_bdt->SetMarkerColor(kMagenta);
        gr_mc_maxrms_and_nbarlayer13_eff_bdt->SetMarkerColor(kGreen);

        gr_data_maxrms_and_nbarlayer13_eff_xtrl_tight->Draw("AP");
        gr_mc_maxrms_and_nbarlayer13_eff_xtrl_tight->Draw("P");
        gr_data_maxrms_and_nbarlayer13_eff_bdt->Draw("P");
        gr_mc_maxrms_and_nbarlayer13_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_maxrms_and_nbarlayer13_eff_xtrl_tight->SetMinimum(0.9);
        gr_data_maxrms_and_nbarlayer13_eff_xtrl_tight->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "maxrms and nbarlayer13 Efficiency (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf","Title:maxrms and nbarlayer13 - xtrl tight efficiency");
        
        // Print maxrms and nbarlayer13 efficiency (xtrl loose)
        auto gr_data_maxrms_and_nbarlayer13_eff_xtrl_loose = remove_x_error_from_TEfficiency(data_maxrms_and_nbarlayer13_eff_xtrl_loose);
        auto gr_mc_maxrms_and_nbarlayer13_eff_xtrl_loose = remove_x_error_from_TEfficiency(mc_maxrms_and_nbarlayer13_eff_xtrl_loose);
        
        gr_data_maxrms_and_nbarlayer13_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_maxrms_and_nbarlayer13_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_maxrms_and_nbarlayer13_eff_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_maxrms_and_nbarlayer13_eff_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_maxrms_and_nbarlayer13_eff_xtrl_loose->Draw("AP");
        gr_mc_maxrms_and_nbarlayer13_eff_xtrl_loose->Draw("P");
        gr_data_maxrms_and_nbarlayer13_eff_bdt->Draw("P");
        gr_mc_maxrms_and_nbarlayer13_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_maxrms_and_nbarlayer13_eff_xtrl_loose->SetMinimum(0.9);
        gr_data_maxrms_and_nbarlayer13_eff_xtrl_loose->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "maxrms and nbarlayer13 Efficiency (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:maxrms and nbarlayer13 - xtrl loose efficiency");

        // Print sumrms low energy efficiency (xtrl tight)
        auto gr_data_sumrms_low_energy_eff_xtrl_tight = remove_x_error_from_TEfficiency(data_sumrms_low_energy_eff_xtrl_tight);
        auto gr_mc_sumrms_low_energy_eff_xtrl_tight = remove_x_error_from_TEfficiency(mc_sumrms_low_energy_eff_xtrl_tight);
        auto gr_data_sumrms_low_energy_eff_bdt = remove_x_error_from_TEfficiency(data_sumrms_low_energy_eff_bdt);
        auto gr_mc_sumrms_low_energy_eff_bdt = remove_x_error_from_TEfficiency(mc_sumrms_low_energy_eff_bdt);
        
        gr_data_sumrms_low_energy_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_sumrms_low_energy_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_sumrms_low_energy_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_sumrms_low_energy_eff_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_sumrms_low_energy_eff_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_sumrms_low_energy_eff_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_sumrms_low_energy_eff_bdt->SetMarkerColor(kMagenta);
        gr_mc_sumrms_low_energy_eff_bdt->SetMarkerColor(kGreen);

        gr_data_sumrms_low_energy_eff_xtrl_tight->Draw("AP");
        gr_mc_sumrms_low_energy_eff_xtrl_tight->Draw("P");
        gr_data_sumrms_low_energy_eff_bdt->Draw("P");
        gr_mc_sumrms_low_energy_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_sumrms_low_energy_eff_xtrl_tight->SetMinimum(0.9);
        gr_data_sumrms_low_energy_eff_xtrl_tight->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "sumrms low energy Efficiency (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf","Title:sumrms low energy efficiency - xtrl tight efficiency");
        
        // Print sumrms low energy efficiency (xtrl loose)
        auto gr_data_sumrms_low_energy_eff_xtrl_loose = remove_x_error_from_TEfficiency(data_sumrms_low_energy_eff_xtrl_loose);
        auto gr_mc_sumrms_low_energy_eff_xtrl_loose = remove_x_error_from_TEfficiency(mc_sumrms_low_energy_eff_xtrl_loose);
        
        gr_data_sumrms_low_energy_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_sumrms_low_energy_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_sumrms_low_energy_eff_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_sumrms_low_energy_eff_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_sumrms_low_energy_eff_xtrl_tight->Draw("AP");
        gr_mc_sumrms_low_energy_eff_xtrl_tight->Draw("P");
        gr_data_sumrms_low_energy_eff_bdt->Draw("P");
        gr_mc_sumrms_low_energy_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_sumrms_low_energy_eff_xtrl_tight->SetMinimum(0.9);
        gr_data_sumrms_low_energy_eff_xtrl_tight->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "sumrms low energy Efficiency (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:sumrms low energy - xtrl loose efficiency");

        // Print trackselection efficiency (xtrl tight)
        auto gr_data_track_selection_eff_xtrl_tight = remove_x_error_from_TEfficiency(data_track_selection_eff_xtrl_tight);
        auto gr_mc_track_selection_eff_xtrl_tight = remove_x_error_from_TEfficiency(mc_track_selection_eff_xtrl_tight);
        auto gr_data_track_selection_eff_bdt = remove_x_error_from_TEfficiency(data_track_selection_eff_bdt);
        auto gr_mc_track_selection_eff_bdt = remove_x_error_from_TEfficiency(mc_track_selection_eff_bdt);
        
        gr_data_track_selection_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_track_selection_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_track_selection_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_track_selection_eff_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_track_selection_eff_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_track_selection_eff_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_track_selection_eff_bdt->SetMarkerColor(kMagenta);
        gr_mc_track_selection_eff_bdt->SetMarkerColor(kGreen);

        gr_data_track_selection_eff_xtrl_tight->Draw("AP");
        gr_mc_track_selection_eff_xtrl_tight->Draw("P");
        gr_data_track_selection_eff_bdt->Draw("P");
        gr_mc_track_selection_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_track_selection_eff_xtrl_tight->SetMinimum(0.8);
        gr_data_track_selection_eff_xtrl_tight->SetMaximum(1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "trackselection Efficiency (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf","Title:trackselection - xtrl tight efficiency");
        
        // Print trackselection efficiency (xtrl loose)
        auto gr_data_track_selection_eff_xtrl_loose = remove_x_error_from_TEfficiency(data_track_selection_eff_xtrl_loose);
        auto gr_mc_track_selection_eff_xtrl_loose = remove_x_error_from_TEfficiency(mc_track_selection_eff_xtrl_loose);
        
        gr_data_track_selection_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_track_selection_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_track_selection_eff_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_track_selection_eff_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_track_selection_eff_xtrl_loose->Draw("AP");
        gr_mc_track_selection_eff_xtrl_loose->Draw("P");
        gr_data_track_selection_eff_bdt->Draw("P");
        gr_mc_track_selection_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_track_selection_eff_xtrl_loose->SetMinimum(0.8);
        gr_data_track_selection_eff_xtrl_loose->SetMaximum(1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "trackselection Efficiency (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:trackselection - xtrl loose efficiency");

        // Print trackselection efficiency within STK fiducial volume (xtrl tight)
        auto gr_data_track_selection_eff_within_stk_fvolume_xtrl_tight = remove_x_error_from_TEfficiency(data_track_selection_eff_within_stk_fvolume_xtrl_tight);
        auto gr_mc_track_selection_eff_within_stk_fvolume_xtrl_tight = remove_x_error_from_TEfficiency(mc_track_selection_eff_within_stk_fvolume_xtrl_tight);
        auto gr_data_track_selection_eff_within_stk_fvolume_bdt = remove_x_error_from_TEfficiency(data_track_selection_eff_within_stk_fvolume_bdt);
        auto gr_mc_track_selection_eff_within_stk_fvolume_bdt = remove_x_error_from_TEfficiency(mc_track_selection_eff_within_stk_fvolume_bdt);
        
        gr_data_track_selection_eff_within_stk_fvolume_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_track_selection_eff_within_stk_fvolume_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_track_selection_eff_within_stk_fvolume_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_track_selection_eff_within_stk_fvolume_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_track_selection_eff_within_stk_fvolume_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_track_selection_eff_within_stk_fvolume_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_track_selection_eff_within_stk_fvolume_bdt->SetMarkerColor(kMagenta);
        gr_mc_track_selection_eff_within_stk_fvolume_bdt->SetMarkerColor(kGreen);

        gr_data_track_selection_eff_within_stk_fvolume_xtrl_tight->Draw("AP");
        gr_mc_track_selection_eff_within_stk_fvolume_xtrl_tight->Draw("P");
        gr_data_track_selection_eff_within_stk_fvolume_bdt->Draw("P");
        gr_mc_track_selection_eff_within_stk_fvolume_bdt->Draw("P");

        gPad->Update(); 
        gr_data_track_selection_eff_within_stk_fvolume_xtrl_tight->SetMinimum(0.9);
        gr_data_track_selection_eff_within_stk_fvolume_xtrl_tight->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "trackselection Efficiency within STK fiducial volume (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf","Title:trackselection within STK fiducial volume - xtrl tight efficiency");
        
        // Print trackselection efficiency within STK fiducial volume (xtrl loose)
        auto gr_data_track_selection_eff_within_stk_fvolume_xtrl_loose = remove_x_error_from_TEfficiency(data_track_selection_eff_within_stk_fvolume_xtrl_loose);
        auto gr_mc_track_selection_eff_within_stk_fvolume_xtrl_loose = remove_x_error_from_TEfficiency(mc_track_selection_eff_within_stk_fvolume_xtrl_loose);
        
        gr_data_track_selection_eff_within_stk_fvolume_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_track_selection_eff_within_stk_fvolume_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_track_selection_eff_within_stk_fvolume_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_track_selection_eff_within_stk_fvolume_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_track_selection_eff_within_stk_fvolume_xtrl_loose->Draw("AP");
        gr_mc_track_selection_eff_within_stk_fvolume_xtrl_loose->Draw("P");
        gr_data_track_selection_eff_within_stk_fvolume_bdt->Draw("P");
        gr_mc_track_selection_eff_within_stk_fvolume_bdt->Draw("P");

        gPad->Update(); 
        gr_data_track_selection_eff_within_stk_fvolume_xtrl_loose->SetMinimum(0.9);
        gr_data_track_selection_eff_within_stk_fvolume_xtrl_loose->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "trackselection Efficiency within STK fiducial volume (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:trackselection within STK fiducial volume - xtrl loose efficiency");

        // Print STK 1 RM efficiency (xtrl tight)
        auto gr_data_stk_1_rm_eff_xtrl_tight = remove_x_error_from_TEfficiency(data_stk_1_rm_eff_xtrl_tight);
        auto gr_mc_stk_1_rm_eff_xtrl_tight = remove_x_error_from_TEfficiency(mc_stk_1_rm_eff_xtrl_tight);
        auto gr_data_stk_1_rm_eff_bdt = remove_x_error_from_TEfficiency(data_stk_1_rm_eff_bdt);
        auto gr_mc_stk_1_rm_eff_bdt = remove_x_error_from_TEfficiency(mc_stk_1_rm_eff_bdt);
        
        gr_data_stk_1_rm_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_stk_1_rm_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_stk_1_rm_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_stk_1_rm_eff_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_stk_1_rm_eff_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_stk_1_rm_eff_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_stk_1_rm_eff_bdt->SetMarkerColor(kMagenta);
        gr_mc_stk_1_rm_eff_bdt->SetMarkerColor(kGreen);

        gr_data_stk_1_rm_eff_xtrl_tight->Draw("AP");
        gr_mc_stk_1_rm_eff_xtrl_tight->Draw("P");
        gr_data_stk_1_rm_eff_bdt->Draw("P");
        gr_mc_stk_1_rm_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_stk_1_rm_eff_xtrl_tight->SetMinimum(0.9);
        gr_data_stk_1_rm_eff_xtrl_tight->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK 1 RM efficiency (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf","Title:STK 1 RM efficiency - xtrl tight efficiency");
        
        // Print STK 1 RM efficiency (xtrl loose)
        auto gr_data_stk_1_rm_eff_xtrl_loose = remove_x_error_from_TEfficiency(data_stk_1_rm_eff_xtrl_loose);
        auto gr_mc_stk_1_rm_eff_xtrl_loose = remove_x_error_from_TEfficiency(mc_stk_1_rm_eff_xtrl_loose);
        
        gr_data_stk_1_rm_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_stk_1_rm_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_stk_1_rm_eff_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_stk_1_rm_eff_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_stk_1_rm_eff_xtrl_loose->Draw("AP");
        gr_mc_stk_1_rm_eff_xtrl_loose->Draw("P");
        gr_data_stk_1_rm_eff_bdt->Draw("P");
        gr_mc_stk_1_rm_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_stk_1_rm_eff_xtrl_loose->SetMinimum(0.9);
        gr_data_stk_1_rm_eff_xtrl_loose->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK 1 RM efficiency (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:STK 1 RM efficiency - xtrl loose efficiency");

        // Print psd-stk match efficiency (xtrl tight)
        auto gr_data_psd_stk_match_eff_xtrl_tight = remove_x_error_from_TEfficiency(data_psd_stk_match_eff_xtrl_tight);
        auto gr_mc_psd_stk_match_eff_xtrl_tight = remove_x_error_from_TEfficiency(mc_psd_stk_match_eff_xtrl_tight);
        auto gr_data_psd_stk_match_eff_bdt = remove_x_error_from_TEfficiency(data_psd_stk_match_eff_bdt);
        auto gr_mc_psd_stk_match_eff_bdt = remove_x_error_from_TEfficiency(mc_psd_stk_match_eff_bdt);

        gr_data_psd_stk_match_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_psd_stk_match_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_psd_stk_match_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_psd_stk_match_eff_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_psd_stk_match_eff_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_psd_stk_match_eff_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_psd_stk_match_eff_bdt->SetMarkerColor(kMagenta);
        gr_mc_psd_stk_match_eff_bdt->SetMarkerColor(kGreen);

        gr_data_psd_stk_match_eff_xtrl_tight->Draw("AP");
        gr_mc_psd_stk_match_eff_xtrl_tight->Draw("P");
        gr_data_psd_stk_match_eff_bdt->Draw("P");
        gr_mc_psd_stk_match_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_psd_stk_match_eff_xtrl_tight->SetMinimum(0.9);
        gr_data_psd_stk_match_eff_xtrl_tight->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "psd-stk match Efficiency (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf","Title:psd-stk match - xtrl tight efficiency");
        
        // Print psd-stk match efficiency (xtrl loose)
        auto gr_data_psd_stk_match_eff_xtrl_loose = remove_x_error_from_TEfficiency(data_psd_stk_match_eff_xtrl_loose);
        auto gr_mc_psd_stk_match_eff_xtrl_loose = remove_x_error_from_TEfficiency(mc_psd_stk_match_eff_xtrl_loose);
        
        gr_data_psd_stk_match_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_psd_stk_match_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_psd_stk_match_eff_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_psd_stk_match_eff_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_psd_stk_match_eff_xtrl_loose->Draw("AP");
        gr_mc_psd_stk_match_eff_xtrl_loose->Draw("P");
        gr_data_psd_stk_match_eff_bdt->Draw("P");
        gr_mc_psd_stk_match_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_psd_stk_match_eff_xtrl_loose->SetMinimum(0.9);
        gr_data_psd_stk_match_eff_xtrl_loose->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "psd-stk match Efficiency (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:psd-stk match - xtrl loose efficiency");

        // Print psd charge efficiency (xtrl tight)
        auto gr_data_psd_charge_eff_xtrl_tight = remove_x_error_from_TEfficiency(data_psd_charge_eff_xtrl_tight);
        auto gr_mc_psd_charge_eff_xtrl_tight = remove_x_error_from_TEfficiency(mc_psd_charge_eff_xtrl_tight);
        auto gr_data_psd_charge_eff_bdt = remove_x_error_from_TEfficiency(data_psd_charge_eff_bdt);
        auto gr_mc_psd_charge_eff_bdt = remove_x_error_from_TEfficiency(mc_psd_charge_eff_bdt);
        
        gr_data_psd_charge_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_psd_charge_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_psd_charge_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_psd_charge_eff_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_psd_charge_eff_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_psd_charge_eff_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_psd_charge_eff_bdt->SetMarkerColor(kMagenta);
        gr_mc_psd_charge_eff_bdt->SetMarkerColor(kGreen);

        gr_data_psd_charge_eff_xtrl_tight->Draw("AP");
        gr_mc_psd_charge_eff_xtrl_tight->Draw("P");
        gr_data_psd_charge_eff_bdt->Draw("P");
        gr_mc_psd_charge_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_psd_charge_eff_xtrl_tight->SetMinimum(0.9);
        gr_data_psd_charge_eff_xtrl_tight->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "PSD charge Efficiency (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf","Title:PSD charge - xtrl tight efficiency");
        
        // Print psd charge efficiency (xtrl loose)
        auto gr_data_psd_charge_eff_xtrl_loose = remove_x_error_from_TEfficiency(data_psd_charge_eff_xtrl_loose);
        auto gr_mc_psd_charge_eff_xtrl_loose = remove_x_error_from_TEfficiency(mc_psd_charge_eff_xtrl_loose);
        
        gr_data_psd_charge_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_psd_charge_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_psd_charge_eff_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_psd_charge_eff_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_psd_charge_eff_xtrl_loose->Draw("AP");
        gr_mc_psd_charge_eff_xtrl_loose->Draw("P");
        gr_data_psd_charge_eff_bdt->Draw("P");
        gr_mc_psd_charge_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_psd_charge_eff_xtrl_loose->SetMinimum(0.9);
        gr_data_psd_charge_eff_xtrl_loose->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "PSD charge Efficiency (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:PSD charge - xtrl loose efficiency");

        // Print stk charge efficiency (xtrl tight)
        auto gr_data_stk_charge_eff_xtrl_tight = remove_x_error_from_TEfficiency(data_stk_charge_eff_xtrl_tight);
        auto gr_mc_stk_charge_eff_xtrl_tight = remove_x_error_from_TEfficiency(mc_stk_charge_eff_xtrl_tight);
        auto gr_data_stk_charge_eff_bdt = remove_x_error_from_TEfficiency(data_stk_charge_eff_bdt);
        auto gr_mc_stk_charge_eff_bdt = remove_x_error_from_TEfficiency(mc_stk_charge_eff_bdt);
        
        gr_data_stk_charge_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_stk_charge_eff_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_stk_charge_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_stk_charge_eff_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_stk_charge_eff_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_stk_charge_eff_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_stk_charge_eff_bdt->SetMarkerColor(kMagenta);
        gr_mc_stk_charge_eff_bdt->SetMarkerColor(kGreen);

        gr_data_stk_charge_eff_xtrl_tight->Draw("AP");
        gr_mc_stk_charge_eff_xtrl_tight->Draw("P");
        gr_data_stk_charge_eff_bdt->Draw("P");
        gr_mc_stk_charge_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_stk_charge_eff_xtrl_tight->SetMinimum(0.95);
        gr_data_stk_charge_eff_xtrl_tight->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK charge Efficiency (xtrl tight)", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency.pdf","Title:STK charge - xtrl tight efficiency");
        
        // Print stk charge efficiency (xtrl loose)
        auto gr_data_stk_charge_eff_xtrl_loose = remove_x_error_from_TEfficiency(data_stk_charge_eff_xtrl_loose);
        auto gr_mc_stk_charge_eff_xtrl_loose = remove_x_error_from_TEfficiency(mc_stk_charge_eff_xtrl_loose);
        
        gr_data_stk_charge_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_stk_charge_eff_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        
        gr_data_stk_charge_eff_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_stk_charge_eff_xtrl_loose->SetMarkerColor(kBlue);
        
        gr_data_stk_charge_eff_xtrl_loose->Draw("AP");
        gr_mc_stk_charge_eff_xtrl_loose->Draw("P");
        gr_data_stk_charge_eff_bdt->Draw("P");
        gr_mc_stk_charge_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_stk_charge_eff_xtrl_loose->SetMinimum(0.95);
        gr_data_stk_charge_eff_xtrl_loose->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK charge Efficiency (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:STK charge - xtrl loose efficiency");

        // Print STK cluster efficiency on first layer within PSD fvolume (xtrl tight)
        auto gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight = remove_x_error_from_TEfficiency(data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight);
        auto gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight = remove_x_error_from_TEfficiency(mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight);
        auto gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt = remove_x_error_from_TEfficiency(data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt);
        auto gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt = remove_x_error_from_TEfficiency(mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt);
        
        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt->SetMarkerColor(kMagenta);
        gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt->SetMarkerColor(kGreen);

        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight->Draw("AP");
        gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight->Draw("P");
        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt->Draw("P");
        gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt->Draw("P");

        gPad->Update(); 
        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight->SetMinimum(0);
        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_tight->SetMaximum(1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK cluster efficiency on first layer within PSD fvolume (xtrl tight)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:STK cluster efficiency on first layer within PSD fvolume");

        // Print STK cluster efficiency on first layer within PSD fvolume (xtrl loose)
        auto gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose = remove_x_error_from_TEfficiency(data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose);
        auto gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose = remove_x_error_from_TEfficiency(mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose);
        
        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose->Draw("AP");
        gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose->Draw("P");
        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt->Draw("P");
        gr_mc_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_bdt->Draw("P");

        gPad->Update(); 
        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose->SetMinimum(0);
        gr_data_clusters_on_first_stk_layer_eff_within_psd_fvolume_xtrl_loose->SetMaximum(1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK cluster efficiency on first layer within PSD fvolume (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:STK cluster efficiency on first layer within PSD fvolume");

        // Print STK cluster efficiency on first layer outside PSD fvolume (xtrl tight)
        auto gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight = remove_x_error_from_TEfficiency(data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight);
        auto gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight = remove_x_error_from_TEfficiency(mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight);
        auto gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt = remove_x_error_from_TEfficiency(data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt);
        auto gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt = remove_x_error_from_TEfficiency(mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt);
        
        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight->SetMarkerStyle(kFullDotMedium);
        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight->SetMarkerColor(kRed);
        gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight->SetMarkerColor(kBlue);
        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt->SetMarkerColor(kMagenta);
        gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt->SetMarkerColor(kGreen);

        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight->Draw("AP");
        gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight->Draw("P");
        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt->Draw("P");
        gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt->Draw("P");

        gPad->Update(); 
        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight->SetMinimum(0);
        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_tight->SetMaximum(1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK cluster efficiency on first layer outside PSD fvolume (xtrl tight)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf","Title:STK cluster efficiency on first layer outside PSD fvolume");

        // Print STK cluster efficiency on first layer outside PSD fvolume (xtrl loose)
        auto gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose = remove_x_error_from_TEfficiency(data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose);
        auto gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose = remove_x_error_from_TEfficiency(mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose);

        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose->SetMarkerStyle(kFullDotMedium);
        gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose->SetMarkerStyle(kFullDotMedium);

        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose->SetMarkerColor(kRed);
        gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose->SetMarkerColor(kBlue);

        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose->Draw("AP");
        gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose->Draw("P");
        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt->Draw("P");
        gr_mc_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_bdt->Draw("P");

        gPad->Update(); 
        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose->SetMinimum(0);
        gr_data_clusters_on_first_stk_layer_eff_outside_psd_fvolume_xtrl_loose->SetMaximum(1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK cluster efficiency on first layer outside PSD fvolume (xtrl loose)", "tlNDC");
        label.Draw();
        print_canvas.Print("efficiency.pdf)","Title:STK cluster efficiency on first layer outside PSD fvolume");
        
    }