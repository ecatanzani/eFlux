#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iostream>

#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"

const double generation_vertex_radius = 1.381976597885342;

std::vector<double> getEnergyBinning(const char* config_file);
std::vector<double> createLogBinning(const double eMin, const double eMax, const std::size_t n_bins);
double wtsydp(const double minene, const double maxene, const double index);

void buildAcceptance(
    const char* acc_config_file,
    const char* full_acceptance_histos, 
    const char* output_file = "acceptance.root")
{
    auto energy_binning = getEnergyBinning(acc_config_file);

    TFile* input_file = TFile::Open(full_acceptance_histos, "READ");
    if (input_file->IsZombie())
    {
        std::cerr << "\n\nError opening input file: [" << full_acceptance_histos << "]\n\n";
        exit(100);
    }

    TH1D* h_gen                                     = static_cast<TH1D*>(input_file->Get("h_gen"));
    TH1D* h_geometric                               = static_cast<TH1D*>(input_file->Get("h_geometric"));
    TH1D* h_geometric_trigger                       = static_cast<TH1D*>(input_file->Get("h_geometric_trigger"));
    TH1D* h_bgo_fiducial                            = static_cast<TH1D*>(input_file->Get("h_bgo_fiducial"));
    TH1D* h_bgo_fiducial_het                        = static_cast<TH1D*>(input_file->Get("h_bgo_fiducial_het"));
    TH1D* h_nBarLayer13                             = static_cast<TH1D*>(input_file->Get("h_nBarLayer13"));
    TH1D* h_maxrms                                  = static_cast<TH1D*>(input_file->Get("h_maxrms"));
    TH1D* h_trackselection                          = static_cast<TH1D*>(input_file->Get("h_trackselection"));
    TH1D* h_trackselection_no_3hit_recover          = static_cast<TH1D*>(input_file->Get("h_trackselection_no_3hit_recover"));
    TH1D* h_bgo_stk_selection_inside_psd_fvolume    = static_cast<TH1D*>(input_file->Get("h_bgo_stk_selection_inside_psd_fvolume"));
    TH1D* h_bgo_stk_selection_outside_psd_fvolume   = static_cast<TH1D*>(input_file->Get("h_bgo_stk_selection_outside_psd_fvolume"));
    TH1D* h_psdstkmatch                             = static_cast<TH1D*>(input_file->Get("h_psdstkmatch"));
    TH1D* h_psdcharge_no_one_view_recover           = static_cast<TH1D*>(input_file->Get("h_psdcharge_no_one_view_recover"));
    TH1D* h_all_cut                                 = static_cast<TH1D*>(input_file->Get("h_all_cut"));

    h_gen                                       ->SetDirectory(0);
    h_geometric                                 ->SetDirectory(0);
    h_geometric_trigger                         ->SetDirectory(0);
    h_bgo_fiducial                              ->SetDirectory(0);
    h_bgo_fiducial_het                          ->SetDirectory(0);
    h_nBarLayer13                               ->SetDirectory(0);
    h_maxrms                                    ->SetDirectory(0);
    h_trackselection                            ->SetDirectory(0);
    h_trackselection_no_3hit_recover            ->SetDirectory(0);
    h_bgo_stk_selection_inside_psd_fvolume      ->SetDirectory(0);
    h_bgo_stk_selection_outside_psd_fvolume     ->SetDirectory(0);
    h_psdstkmatch                               ->SetDirectory(0);
    h_psdcharge_no_one_view_recover             ->SetDirectory(0);
    h_all_cut                                   ->SetDirectory(0);

    input_file->Close();

    TH1D* h_geometric_factor                            = static_cast<TH1D*>(h_geometric->Clone("h_geometric_factor"));
    TH1D* h_acc_geometric_trigger                       = static_cast<TH1D*>(h_geometric_trigger->Clone("h_acc_geometric_trigger"));
    TH1D* h_acc_bgo_fiducial                            = static_cast<TH1D*>(h_bgo_fiducial->Clone("h_acc_bgo_fiducial"));
    TH1D* h_acc_bgo_fiducial_het                        = static_cast<TH1D*>(h_bgo_fiducial_het->Clone("h_acc_bgo_fiducial_het"));
    TH1D* h_acc_nBarLayer13                             = static_cast<TH1D*>(h_nBarLayer13->Clone("h_acc_nBarLayer13"));
    TH1D* h_acc_maxrms                                  = static_cast<TH1D*>(h_maxrms->Clone("h_acc_maxrms"));
    TH1D* h_acc_trackselection                          = static_cast<TH1D*>(h_trackselection->Clone("h_acc_trackselection"));
    TH1D* h_acc_trackselection_no_3hit_recover          = static_cast<TH1D*>(h_trackselection_no_3hit_recover->Clone("h_acc_trackselection_no_3hit_recover"));
    TH1D* h_acc_bgo_stk_selection_inside_psd_fvolume    = static_cast<TH1D*>(h_bgo_stk_selection_inside_psd_fvolume->Clone("h_acc_bgo_stk_selection_inside_psd_fvolume"));
    TH1D* h_acc_bgo_stk_selection_outside_psd_fvolume   = static_cast<TH1D*>(h_bgo_stk_selection_outside_psd_fvolume->Clone("h_acc_bgo_stk_selection_outside_psd_fvolume"));
    TH1D* h_acc_psdstkmatch                             = static_cast<TH1D*>(h_psdstkmatch->Clone("h_acc_psdstkmatch"));
    TH1D* h_acc_psdcharge_no_one_view_recover           = static_cast<TH1D*>(h_psdcharge_no_one_view_recover->Clone("h_acc_psdcharge_no_one_view_recover"));
    TH1D* h_acc_all_cut                                 = static_cast<TH1D*>(h_all_cut->Clone("h_acc_all_cut"));

    h_geometric_factor                              ->Divide(h_gen);
    h_acc_geometric_trigger                         ->Divide(h_gen);
    h_acc_bgo_fiducial                              ->Divide(h_gen);
    h_acc_bgo_fiducial_het                          ->Divide(h_gen);
    h_acc_nBarLayer13                               ->Divide(h_gen);
    h_acc_maxrms                                    ->Divide(h_gen);
    h_acc_trackselection                            ->Divide(h_gen);
    h_acc_trackselection_no_3hit_recover            ->Divide(h_gen);
    h_acc_bgo_stk_selection_inside_psd_fvolume      ->Divide(h_gen);
    h_acc_bgo_stk_selection_outside_psd_fvolume     ->Divide(h_gen);
    h_acc_psdstkmatch                               ->Divide(h_gen);
    h_acc_psdcharge_no_one_view_recover             ->Divide(h_gen);
    h_acc_all_cut                                   ->Divide(h_gen);

    double genSurface               = 4 * TMath::Pi() * pow(generation_vertex_radius, 2) / 2;
    double scaleFactor              = TMath::Pi() * genSurface;

    h_geometric_factor                              ->Scale(scaleFactor);
    h_acc_geometric_trigger                         ->Scale(scaleFactor);
    h_acc_bgo_fiducial                              ->Scale(scaleFactor);
    h_acc_bgo_fiducial_het                          ->Scale(scaleFactor);
    h_acc_nBarLayer13                               ->Scale(scaleFactor);
    h_acc_maxrms                                    ->Scale(scaleFactor);
    h_acc_trackselection                            ->Scale(scaleFactor);
    h_acc_trackselection_no_3hit_recover            ->Scale(scaleFactor);
    h_acc_bgo_stk_selection_inside_psd_fvolume      ->Scale(scaleFactor);
    h_acc_bgo_stk_selection_outside_psd_fvolume     ->Scale(scaleFactor);
    h_acc_psdstkmatch                               ->Scale(scaleFactor);
    h_acc_psdcharge_no_one_view_recover             ->Scale(scaleFactor);
    h_acc_all_cut                                   ->Scale(scaleFactor);

    h_geometric_factor                              ->GetYaxis()->SetTitle("geometric factor [m^{2} sr]");
    h_acc_geometric_trigger                         ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_bgo_fiducial                              ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_bgo_fiducial_het                          ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_nBarLayer13                               ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_maxrms                                    ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_trackselection                            ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_trackselection_no_3hit_recover            ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_bgo_stk_selection_inside_psd_fvolume      ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_bgo_stk_selection_outside_psd_fvolume     ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_psdstkmatch                               ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_psdcharge_no_one_view_recover             ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_all_cut                                   ->GetYaxis()->SetTitle("acceptance [m^{2} sr]");

    std::vector<double> geometric_factor                            (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_geometric_trigger                       (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_fiducial                            (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_fiducial_het                        (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_nBarLayer13                             (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_maxrms                                  (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_trackselection                          (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_trackselection_no_3hit_recover          (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_stk_selection_inside_psd_fvolume    (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_stk_selection_outside_psd_fvolume   (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_psdstkmatch                             (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_psdcharge_no_one_view_recover           (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_all_cut                                 (h_geometric_factor->GetNbinsX(), 0);

    std::vector<double> geometric_factor_err                            (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_geometric_trigger_err                       (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_fiducial_err                            (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_fiducial_het_err                        (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_nBarLayer13_err                             (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_maxrms_err                                  (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_trackselection_err                          (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_trackselection_no_3hit_recover_err          (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_stk_selection_inside_psd_fvolume_err    (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_stk_selection_outside_psd_fvolume_err   (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_psdstkmatch_err                             (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_psdcharge_no_one_view_recover_err           (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_all_cut_err                                 (h_geometric_factor->GetNbinsX(), 0);

    std::vector<double> energy                                          (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> energy_err                                      (h_geometric_factor->GetNbinsX(), 0);

    std::vector<double> ratio_trackselection_no_3hit_recover_trackselection     (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> ratio_psdcharge_no_one_view_recover_psdcharge           (h_geometric_factor->GetNbinsX(), 0);

    for (unsigned int idx=0; idx<energy.size(); ++idx)
    {
        energy[idx]                                             = wtsydp(energy_binning[idx], energy_binning[idx+1], -2.7);
        geometric_factor[idx]                                   = h_geometric_factor->GetBinContent(idx+1);
        geometric_factor_err[idx]                               = h_geometric_factor->GetBinError(idx+1);
        acc_geometric_trigger[idx]                              = h_acc_geometric_trigger->GetBinContent(idx+1);
        acc_geometric_trigger_err[idx]                          = h_acc_geometric_trigger->GetBinError(idx+1);
        acc_bgo_fiducial[idx]                                   = h_acc_bgo_fiducial->GetBinContent(idx+1);
        acc_bgo_fiducial_err[idx]                               = h_acc_bgo_fiducial->GetBinError(idx+1);
        acc_bgo_fiducial_het[idx]                               = h_acc_bgo_fiducial_het->GetBinContent(idx+1);
        acc_bgo_fiducial_het_err[idx]                           = h_acc_bgo_fiducial_het->GetBinError(idx+1);
        acc_nBarLayer13[idx]                                    = h_acc_nBarLayer13->GetBinContent(idx+1);
        acc_nBarLayer13_err[idx]                                = h_acc_nBarLayer13->GetBinError(idx+1);
        acc_maxrms[idx]                                         = h_acc_maxrms->GetBinContent(idx+1);
        acc_maxrms_err[idx]                                     = h_acc_maxrms->GetBinError(idx+1);
        acc_trackselection[idx]                                 = h_acc_trackselection->GetBinContent(idx+1);
        acc_trackselection_err[idx]                             = h_acc_trackselection->GetBinError(idx+1);
        acc_trackselection_no_3hit_recover[idx]                 = h_acc_trackselection_no_3hit_recover->GetBinContent(idx+1);
        acc_trackselection_no_3hit_recover_err[idx]             = h_acc_trackselection_no_3hit_recover->GetBinError(idx+1);
        acc_bgo_stk_selection_inside_psd_fvolume[idx]           = h_acc_bgo_stk_selection_inside_psd_fvolume->GetBinContent(idx+1);
        acc_bgo_stk_selection_inside_psd_fvolume_err[idx]       = h_acc_bgo_stk_selection_inside_psd_fvolume->GetBinError(idx+1);
        acc_bgo_stk_selection_outside_psd_fvolume[idx]          = h_acc_bgo_stk_selection_outside_psd_fvolume->GetBinContent(idx+1);
        acc_bgo_stk_selection_outside_psd_fvolume_err[idx]      = h_acc_bgo_stk_selection_outside_psd_fvolume->GetBinError(idx+1);
        acc_psdstkmatch[idx]                                    = h_acc_psdstkmatch->GetBinContent(idx+1);
        acc_psdstkmatch_err[idx]                                = h_acc_psdstkmatch->GetBinError(idx+1);
        acc_psdcharge_no_one_view_recover[idx]                  = h_acc_psdcharge_no_one_view_recover->GetBinContent(idx+1);
        acc_psdcharge_no_one_view_recover_err[idx]              = h_acc_psdcharge_no_one_view_recover->GetBinError(idx+1);
        acc_all_cut[idx]                                        = h_acc_all_cut->GetBinContent(idx+1);
        acc_all_cut_err[idx]                                    = h_acc_all_cut->GetBinError(idx+1);

        ratio_trackselection_no_3hit_recover_trackselection[idx] = acc_trackselection[idx] ? (acc_trackselection_no_3hit_recover[idx]/acc_trackselection[idx])*100 : 0;
        ratio_psdcharge_no_one_view_recover_psdcharge[idx]       = acc_all_cut[idx] ? (acc_psdcharge_no_one_view_recover[idx]/acc_all_cut[idx])*100 : 0;
    }

    TGraphErrors gr_geometric_factor                            ((int)geometric_factor.size(), &energy[0], &geometric_factor[0], &energy_err[0], &geometric_factor_err[0]);
    TGraphErrors gr_acc_geometric_trigger                       ((int)geometric_factor.size(), &energy[0], &acc_geometric_trigger[0], &energy_err[0], &acc_geometric_trigger_err[0]);
    TGraphErrors gr_acc_bgo_fiducial                            ((int)geometric_factor.size(), &energy[0], &acc_bgo_fiducial[0], &energy_err[0], &acc_bgo_fiducial_err[0]);
    TGraphErrors gr_acc_bgo_fiducial_het                        ((int)geometric_factor.size(), &energy[0], &acc_bgo_fiducial_het[0], &energy_err[0], &acc_bgo_fiducial_het_err[0]);
    TGraphErrors gr_acc_nBarLayer13                             ((int)geometric_factor.size(), &energy[0], &acc_nBarLayer13[0], &energy_err[0], &acc_nBarLayer13_err[0]);
    TGraphErrors gr_acc_maxrms                                  ((int)geometric_factor.size(), &energy[0], &acc_maxrms[0], &energy_err[0], &acc_maxrms_err[0]);
    TGraphErrors gr_acc_trackselection                          ((int)geometric_factor.size(), &energy[0], &acc_trackselection[0], &energy_err[0], &acc_trackselection_err[0]);
    TGraphErrors gr_acc_trackselection_no_3hit_recover          ((int)geometric_factor.size(), &energy[0], &acc_trackselection_no_3hit_recover[0], &energy_err[0], &acc_trackselection_no_3hit_recover_err[0]);
    TGraphErrors gr_acc_bgo_stk_selection_inside_psd_fvolume    ((int)geometric_factor.size(), &energy[0], &acc_bgo_stk_selection_inside_psd_fvolume[0], &energy_err[0], &acc_bgo_stk_selection_inside_psd_fvolume_err[0]);
    TGraphErrors gr_acc_bgo_stk_selection_outside_psd_fvolume   ((int)geometric_factor.size(), &energy[0], &acc_bgo_stk_selection_outside_psd_fvolume[0], &energy_err[0], &acc_bgo_stk_selection_outside_psd_fvolume_err[0]);
    TGraphErrors gr_acc_psdstkmatch                             ((int)geometric_factor.size(), &energy[0], &acc_psdstkmatch[0], &energy_err[0], &acc_psdstkmatch_err[0]);
    TGraphErrors gr_acc_psdcharge_no_one_view_recover           ((int)geometric_factor.size(), &energy[0], &acc_psdcharge_no_one_view_recover[0], &energy_err[0], &acc_psdcharge_no_one_view_recover_err[0]);
    TGraphErrors gr_acc_all_cut                                 ((int)geometric_factor.size(), &energy[0], &acc_all_cut[0], &energy_err[0], &acc_all_cut_err[0]);

    TGraph gr_ratio_trackselection_no_3hit_recover_trackselection (int(energy.size()), &energy[0], &ratio_trackselection_no_3hit_recover_trackselection[0]);
    TGraph gr_ratio_psdcharge_no_one_view_recover_psdcharge       (int(energy.size()), &energy[0], &ratio_psdcharge_no_one_view_recover_psdcharge[0]);

    gr_geometric_factor                             .SetName("gr_geometric_factor");
    gr_acc_geometric_trigger                        .SetName("gr_acc_geometric_trigger");
    gr_acc_bgo_fiducial                             .SetName("gr_acc_bgo_fiducial");
    gr_acc_bgo_fiducial_het                         .SetName("gr_acc_bgo_fiducial_het");
    gr_acc_nBarLayer13                              .SetName("gr_acc_nBarLayer13");
    gr_acc_maxrms                                   .SetName("gr_acc_maxrms");
    gr_acc_trackselection                           .SetName("gr_acc_trackselection");
    gr_acc_trackselection_no_3hit_recover           .SetName("gr_acc_trackselection_no_3hit_recover");
    gr_acc_bgo_stk_selection_inside_psd_fvolume     .SetName("gr_acc_bgo_stk_selection_inside_psd_fvolume");
    gr_acc_bgo_stk_selection_outside_psd_fvolume    .SetName("gr_acc_bgo_stk_selection_outside_psd_fvolume");
    gr_acc_psdstkmatch                              .SetName("gr_acc_psdstkmatch");
    gr_acc_psdcharge_no_one_view_recover            .SetName("gr_acc_psdcharge_no_one_view_recover");
    gr_acc_all_cut                                  .SetName("gr_acc_all_cut");

    gr_ratio_trackselection_no_3hit_recover_trackselection      .SetName("gr_ratio_trackselection_no_3hit_recover_trackselection");
    gr_ratio_psdcharge_no_one_view_recover_psdcharge            .SetName("gr_ratio_psdcharge_no_one_view_recover_psdcharge");

    gr_geometric_factor                             .SetTitle("Geometric Factor");
    gr_acc_geometric_trigger                        .SetTitle("Acceptance - geometric + trigger");
    gr_acc_bgo_fiducial                             .SetTitle("Acceptance - BGO fiducial");
    gr_acc_bgo_fiducial_het                         .SetTitle("Acceptance - BGO fiducial + HET");
    gr_acc_nBarLayer13                              .SetTitle("Acceptance - nBarLayer13");
    gr_acc_maxrms                                   .SetTitle("Acceptance - maxrms");
    gr_acc_trackselection                           .SetTitle("Acceptance - track selection");
    gr_acc_trackselection_no_3hit_recover           .SetTitle("Acceptance - track selection - no 3 hit recover");
    gr_acc_bgo_stk_selection_inside_psd_fvolume     .SetTitle("Acceptance - BGO STK selection inside PSD Fiducial Volume");
    gr_acc_bgo_stk_selection_outside_psd_fvolume    .SetTitle("Acceptance - BGO STK selection outside PSD Fiducial Volume");
    gr_acc_psdstkmatch                              .SetTitle("Acceptance - PSD-STK match");
    gr_acc_psdcharge_no_one_view_recover            .SetTitle("Acceptance - PSD charge - no one view recover");
    gr_acc_all_cut                                  .SetTitle("Acceptance - all cuts");

    gr_ratio_trackselection_no_3hit_recover_trackselection      .SetTitle("Ratio - track selection - no 3 hit recover / track selection");
    gr_ratio_psdcharge_no_one_view_recover_psdcharge            .SetTitle("Ratio - PSD charge - no one view recover / PSD charge");

    gr_geometric_factor                             .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_geometric_trigger                        .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_bgo_fiducial                             .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_bgo_fiducial_het                         .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_nBarLayer13                              .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_maxrms                                   .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_trackselection                           .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_trackselection_no_3hit_recover           .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_bgo_stk_selection_inside_psd_fvolume     .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_bgo_stk_selection_outside_psd_fvolume    .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_psdstkmatch                              .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_psdcharge_no_one_view_recover            .GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_all_cut                                  .GetXaxis()->SetTitle("Energy [GeV]");

    gr_ratio_trackselection_no_3hit_recover_trackselection      .GetXaxis()->SetTitle("Energy [GeV]");
    gr_ratio_psdcharge_no_one_view_recover_psdcharge            .GetXaxis()->SetTitle("Energy [GeV]");

    gr_geometric_factor                             .GetYaxis()->SetTitle("geometric factor [m^{2} sr]");
    gr_acc_geometric_trigger                        .GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_bgo_fiducial                             .GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_bgo_fiducial_het                         .GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_nBarLayer13                              .GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_maxrms                                   .GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_trackselection                           .GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_trackselection_no_3hit_recover           .GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_bgo_stk_selection_inside_psd_fvolume     .GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_bgo_stk_selection_outside_psd_fvolume    .GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_psdstkmatch                              .GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_psdcharge_no_one_view_recover            .GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_all_cut                                  .GetYaxis()->SetTitle("acceptance [m^{2} sr]");

    gr_ratio_trackselection_no_3hit_recover_trackselection      .GetYaxis()->SetTitle("track selection_{no 3 hit recover} / track selection (percentage)");
    gr_ratio_psdcharge_no_one_view_recover_psdcharge            .GetYaxis()->SetTitle("PSD charge_{no one view recover} / PSD charge (percentage)");

    TFile* outfile = TFile::Open(output_file, "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "\n\nError writing output acceptance file: [" << output_file << "]\n\n";
        exit(100);
    }

    h_geometric_factor                              ->Write();
    h_acc_geometric_trigger                         ->Write();
    h_acc_bgo_fiducial                              ->Write();
    h_acc_bgo_fiducial_het                          ->Write();
    h_acc_nBarLayer13                               ->Write();
    h_acc_maxrms                                    ->Write();
    h_acc_trackselection                            ->Write();
    h_acc_trackselection_no_3hit_recover            ->Write();
    h_acc_bgo_stk_selection_inside_psd_fvolume      ->Write();
    h_acc_bgo_stk_selection_outside_psd_fvolume     ->Write();
    h_acc_psdstkmatch                               ->Write();
    h_acc_psdcharge_no_one_view_recover             ->Write();
    h_acc_all_cut                                   ->Write();

    gr_geometric_factor                             .Write();
    gr_acc_geometric_trigger                        .Write();
    gr_acc_bgo_fiducial                             .Write();
    gr_acc_bgo_fiducial_het                         .Write();
    gr_acc_nBarLayer13                              .Write();
    gr_acc_maxrms                                   .Write();
    gr_acc_trackselection                           .Write();
    gr_acc_trackselection_no_3hit_recover           .Write();
    gr_acc_bgo_stk_selection_inside_psd_fvolume     .Write();
    gr_acc_bgo_stk_selection_outside_psd_fvolume    .Write();
    gr_acc_psdstkmatch                              .Write();
    gr_acc_psdcharge_no_one_view_recover            .Write();
    gr_acc_all_cut                                  .Write();

    gr_ratio_trackselection_no_3hit_recover_trackselection      .Write();
    gr_ratio_psdcharge_no_one_view_recover_psdcharge            .Write();

    outfile->Close();

    std::cout << "\n\nOutput file has been written: [" << output_file << "]\n\n";
}

std::vector<double> getEnergyBinning(const char* config_file)
{
    std::size_t n_bins = 0;
    double min_event_energy = -999;
    double max_event_energy = -999;
    std::vector<double> energy_binning;

    std::ifstream input_file(config_file);
    if (!input_file.is_open())
	{
		std::cerr << "\nInput acceptance config file not found [" << config_file << "]\n\n";
		exit(100);
	}
    std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
    std::string tmp_str;
	std::istringstream input_stream(input_string);
	std::string::size_type sz;

    while (input_stream >> tmp_str)
	{
		if (!strcmp(tmp_str.c_str(), "n_energy_bins"))
			input_stream >> n_bins;
		if (!strcmp(tmp_str.c_str(), "min_event_energy"))
		{
			input_stream >> tmp_str;
			min_event_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "max_event_energy"))
		{
			input_stream >> tmp_str;
			max_event_energy = stod(tmp_str, &sz);
		}
	}

    return createLogBinning(min_event_energy, max_event_energy, n_bins);
}   

std::vector<double> createLogBinning(
    const double eMin,
	const double eMax,
	const std::size_t n_bins)
{
    std::vector<double> binning(n_bins+1, 0);
	double log_interval = (log10(eMax) - log10(eMin)) / n_bins;
	for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
		binning[bIdx] = pow(10, log10(eMin) + bIdx * log_interval);

	return binning;
}

double wtsydp(const double minene, const double maxene, const double index)
{
    auto dene = maxene - minene;
    if (index != -1)
        return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * dene)), 1. / index);
    else
        return dene / log(maxene / minene);
}