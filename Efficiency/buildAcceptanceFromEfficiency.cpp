#include <vector>
#include <memory>
#include <iostream>

#include "TFile.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include "TGraphErrors.h"

void buildAcceptanceFromEfficiency(
    const char* input_acceptance_file,
    const char* input_efficiency_file,
    const char* output_file) {

        // Read input acceptance file
        TFile *infile_acc = TFile::Open(input_acceptance_file, "READ");
        if (!infile_acc->IsOpen()) {
            std::cerr << "\n\nError reading input acceptance file [" << input_acceptance_file << "]\n\n";
            exit(100);
        }

        auto gr_acc_bgo_fiducial = static_cast<TGraphErrors*>(infile_acc->Get("gr_acc_bgo_fiducial"));
        std::vector<double> gr_acc_bgo_fiducial_energy (gr_acc_bgo_fiducial->GetN(), 0);
        std::vector<double> gr_acc_bgo_fiducial_value (gr_acc_bgo_fiducial->GetN(), 0);

        for (int idx_p=0; idx_p<gr_acc_bgo_fiducial->GetN(); ++idx_p) {
            gr_acc_bgo_fiducial_energy[idx_p]   = gr_acc_bgo_fiducial->GetPointX(idx_p);
            gr_acc_bgo_fiducial_value[idx_p]    = gr_acc_bgo_fiducial->GetPointY(idx_p);
        }

        infile_acc->Close();

        // Read input efficiency file
        TFile *infile_eff = TFile::Open(input_efficiency_file, "READ");
        if (!infile_eff->IsOpen()) {
            std::cerr << "\n\nError reading input efficiency file [" << input_efficiency_file << "]\n\n";
            exit(100);
        }

        auto trigger_eff_het_over_let_xtrl_tight        = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/trigger_eff_het_over_let_xtrl_tight"));
        auto trigger_eff_het_over_unb_xtrl_tight        = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/trigger_eff_het_over_unb_xtrl_tight"));
        auto maxrms_eff_xtrl_tight                      = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/maxrms_eff_xtrl_tight"));
        auto nbarlayer13_eff_xtrl_tight                 = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/nbarlayer13_eff_xtrl_tight"));
        auto maxrms_and_nbarlayer13_eff_xtrl_tight      = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/maxrms_and_nbarlayer13_eff_xtrl_tight"));
        auto track_selection_eff_xtrl_tight             = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/track_selection_eff_xtrl_tight"));
        auto psd_stk_match_eff_xtrl_tight               = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/psd_stk_match_eff_xtrl_tight"));
        auto psd_charge_eff_xtrl_tight                  = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/psd_charge_eff_xtrl_tight"));

        auto trigger_eff_het_over_let_xtrl_loose        = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/trigger_eff_het_over_let_xtrl_loose"));
        auto trigger_eff_het_over_unb_xtrl_loose        = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/trigger_eff_het_over_unb_xtrl_loose"));
        auto maxrms_eff_xtrl_loose                      = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/maxrms_eff_xtrl_loose"));
        auto nbarlayer13_eff_xtrl_loose                 = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/nbarlayer13_eff_xtrl_loose"));
        auto maxrms_and_nbarlayer13_eff_xtrl_loose      = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/maxrms_and_nbarlayer13_eff_xtrl_loose"));
        auto track_selection_eff_xtrl_loose             = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/track_selection_eff_xtrl_loose"));
        auto psd_stk_match_eff_xtrl_loose               = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/psd_stk_match_eff_xtrl_loose"));
        auto psd_charge_eff_xtrl_loose                  = static_cast<TEfficiency*>(infile_eff->Get("efficiencies/psd_charge_eff_xtrl_loose"));

        std::vector<double> trigger_eff_het_over_let_xtrl_tight_value       (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> trigger_eff_het_over_unb_xtrl_tight_value       (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> maxrms_eff_xtrl_tight_value                     (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> nbarlayer13_eff_xtrl_tight_value                (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> maxrms_and_nbarlayer13_eff_xtrl_tight_value     (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> track_selection_eff_xtrl_tight_value            (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> psd_stk_match_eff_xtrl_tight_value              (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> psd_charge_eff_xtrl_tight_value                 (gr_acc_bgo_fiducial_energy.size(), 0);

        std::vector<double> trigger_eff_het_over_let_xtrl_loose_value       (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> trigger_eff_het_over_unb_xtrl_loose_value       (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> maxrms_eff_xtrl_loose_value                     (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> nbarlayer13_eff_xtrl_loose_value                (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> maxrms_and_nbarlayer13_eff_xtrl_loose_value     (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> track_selection_eff_xtrl_loose_value            (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> psd_stk_match_eff_xtrl_loose_value              (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> psd_charge_eff_xtrl_loose_value                 (gr_acc_bgo_fiducial_energy.size(), 0);

        for (int idx_p=0; idx_p<(int)gr_acc_bgo_fiducial_energy.size(); ++idx_p) {
            trigger_eff_het_over_let_xtrl_tight_value[idx_p] = trigger_eff_het_over_let_xtrl_tight->GetEfficiency(trigger_eff_het_over_let_xtrl_tight->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            trigger_eff_het_over_unb_xtrl_tight_value[idx_p] = trigger_eff_het_over_unb_xtrl_tight->GetEfficiency(trigger_eff_het_over_unb_xtrl_tight->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            maxrms_eff_xtrl_tight_value[idx_p] = maxrms_eff_xtrl_tight->GetEfficiency(maxrms_eff_xtrl_tight->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            nbarlayer13_eff_xtrl_tight_value[idx_p] = nbarlayer13_eff_xtrl_tight->GetEfficiency(nbarlayer13_eff_xtrl_tight->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            maxrms_and_nbarlayer13_eff_xtrl_tight_value[idx_p] = maxrms_and_nbarlayer13_eff_xtrl_tight->GetEfficiency(maxrms_and_nbarlayer13_eff_xtrl_tight->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            track_selection_eff_xtrl_tight_value[idx_p] = track_selection_eff_xtrl_tight->GetEfficiency(track_selection_eff_xtrl_tight->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            psd_stk_match_eff_xtrl_tight_value[idx_p] = psd_stk_match_eff_xtrl_tight->GetEfficiency(psd_stk_match_eff_xtrl_tight->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            psd_charge_eff_xtrl_tight_value[idx_p] = psd_charge_eff_xtrl_tight->GetEfficiency(psd_charge_eff_xtrl_tight->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));

            trigger_eff_het_over_let_xtrl_loose_value[idx_p] = trigger_eff_het_over_let_xtrl_loose->GetEfficiency(trigger_eff_het_over_let_xtrl_loose->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            trigger_eff_het_over_unb_xtrl_loose_value[idx_p] = trigger_eff_het_over_unb_xtrl_loose->GetEfficiency(trigger_eff_het_over_unb_xtrl_loose->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            maxrms_eff_xtrl_loose_value[idx_p] = maxrms_eff_xtrl_loose->GetEfficiency(maxrms_eff_xtrl_loose->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            nbarlayer13_eff_xtrl_loose_value[idx_p] = nbarlayer13_eff_xtrl_loose->GetEfficiency(nbarlayer13_eff_xtrl_loose->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            maxrms_and_nbarlayer13_eff_xtrl_loose_value[idx_p] = maxrms_and_nbarlayer13_eff_xtrl_loose->GetEfficiency(maxrms_and_nbarlayer13_eff_xtrl_loose->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            track_selection_eff_xtrl_loose_value[idx_p] = track_selection_eff_xtrl_loose->GetEfficiency(track_selection_eff_xtrl_loose->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            psd_stk_match_eff_xtrl_loose_value[idx_p] = psd_stk_match_eff_xtrl_loose->GetEfficiency(psd_stk_match_eff_xtrl_loose->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
            psd_charge_eff_xtrl_loose_value[idx_p] = psd_charge_eff_xtrl_loose->GetEfficiency(psd_charge_eff_xtrl_loose->FindFixBin(gr_acc_bgo_fiducial_energy[idx_p]));
        }
        
        infile_eff->Close();

        std::vector<double> acceptance_bgo_fiducial_het_fe_xtrl_tight           (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> acceptance_nbarlayer13_fe_xtrl_tight                (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> accepatance_maxrms_fe_xtrl_tight                    (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> acceptance_maxrms_and_nbarlayer13_fe_xtrl_tight     (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> acceptance_track_selection_fe_xtrl_tight            (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> acceptance_psd_stk_match_fe_xtrl_tight              (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> acceptance_psd_charge_fe_xtrl_tight                 (gr_acc_bgo_fiducial_energy.size(), 0);

        std::vector<double> acceptance_bgo_fiducial_het_fe_xtrl_loose           (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> acceptance_nbarlayer13_fe_xtrl_loose                (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> accepatance_maxrms_fe_xtrl_loose                    (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> acceptance_maxrms_and_nbarlayer13_fe_xtrl_loose     (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> acceptance_track_selection_fe_xtrl_loose            (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> acceptance_psd_stk_match_fe_xtrl_loose              (gr_acc_bgo_fiducial_energy.size(), 0);
        std::vector<double> acceptance_psd_charge_fe_xtrl_loose                 (gr_acc_bgo_fiducial_energy.size(), 0);

        for (int idx_p=0; idx_p<(int)gr_acc_bgo_fiducial_energy.size(); ++idx_p) {
            
            acceptance_bgo_fiducial_het_fe_xtrl_tight[idx_p]        = trigger_eff_het_over_let_xtrl_tight_value[idx_p]*gr_acc_bgo_fiducial_value[idx_p];
            acceptance_nbarlayer13_fe_xtrl_tight[idx_p]             = nbarlayer13_eff_xtrl_tight_value[idx_p]*acceptance_bgo_fiducial_het_fe_xtrl_tight[idx_p];
            accepatance_maxrms_fe_xtrl_tight[idx_p]                 = maxrms_eff_xtrl_tight_value[idx_p]*acceptance_nbarlayer13_fe_xtrl_tight[idx_p];
            acceptance_maxrms_and_nbarlayer13_fe_xtrl_tight[idx_p]  = maxrms_and_nbarlayer13_eff_xtrl_tight_value[idx_p]*acceptance_bgo_fiducial_het_fe_xtrl_tight[idx_p];
            acceptance_track_selection_fe_xtrl_tight[idx_p]         = track_selection_eff_xtrl_tight_value[idx_p]*accepatance_maxrms_fe_xtrl_tight[idx_p];
            acceptance_psd_stk_match_fe_xtrl_tight[idx_p]           = psd_stk_match_eff_xtrl_tight_value[idx_p]*acceptance_track_selection_fe_xtrl_tight[idx_p];
            acceptance_psd_charge_fe_xtrl_tight[idx_p]              = psd_charge_eff_xtrl_tight_value[idx_p]*acceptance_psd_stk_match_fe_xtrl_tight[idx_p];

            acceptance_bgo_fiducial_het_fe_xtrl_loose[idx_p]        = trigger_eff_het_over_let_xtrl_loose_value[idx_p]*gr_acc_bgo_fiducial_value[idx_p];
            acceptance_nbarlayer13_fe_xtrl_loose[idx_p]             = nbarlayer13_eff_xtrl_loose_value[idx_p]*acceptance_bgo_fiducial_het_fe_xtrl_loose[idx_p];
            accepatance_maxrms_fe_xtrl_loose[idx_p]                 = maxrms_eff_xtrl_loose_value[idx_p]*acceptance_nbarlayer13_fe_xtrl_loose[idx_p];
            acceptance_maxrms_and_nbarlayer13_fe_xtrl_loose[idx_p]  = maxrms_and_nbarlayer13_eff_xtrl_loose_value[idx_p]*acceptance_bgo_fiducial_het_fe_xtrl_loose[idx_p];
            acceptance_track_selection_fe_xtrl_loose[idx_p]         = track_selection_eff_xtrl_loose_value[idx_p]*accepatance_maxrms_fe_xtrl_loose[idx_p];
            acceptance_psd_stk_match_fe_xtrl_loose[idx_p]           = psd_stk_match_eff_xtrl_loose_value[idx_p]*acceptance_track_selection_fe_xtrl_loose[idx_p];
            acceptance_psd_charge_fe_xtrl_loose[idx_p]              = psd_charge_eff_xtrl_loose_value[idx_p]*acceptance_psd_stk_match_fe_xtrl_loose[idx_p];
        }

        TGraph gr_acceptance_bgo_fiducial_het_fe_xtrl_tight((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_bgo_fiducial_het_fe_xtrl_tight[0]);
        TGraph gr_acceptance_nbarlayer13_fe_xtrl_tight((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_nbarlayer13_fe_xtrl_tight[0]);
        TGraph gr_accepatance_maxrms_fe_xtrl_tight((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &accepatance_maxrms_fe_xtrl_tight[0]);
        TGraph gr_acceptance_maxrms_and_nbarlayer13_fe_xtrl_tight((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_maxrms_and_nbarlayer13_fe_xtrl_tight[0]);
        TGraph gr_acceptance_track_selection_fe_xtrl_tight((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_track_selection_fe_xtrl_tight[0]);
        TGraph gr_acceptance_psd_stk_match_fe_xtrl_tight((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_psd_stk_match_fe_xtrl_tight[0]);
        TGraph gr_acceptance_psd_charge_fe_xtrl_tight((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_psd_charge_fe_xtrl_tight[0]);

        TGraph gr_acceptance_bgo_fiducial_het_fe_xtrl_loose((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_bgo_fiducial_het_fe_xtrl_loose[0]);
        TGraph gr_acceptance_nbarlayer13_fe_xtrl_loose((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_nbarlayer13_fe_xtrl_loose[0]);
        TGraph gr_accepatance_maxrms_fe_xtrl_loose((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &accepatance_maxrms_fe_xtrl_loose[0]);
        TGraph gr_acceptance_maxrms_and_nbarlayer13_fe_xtrl_loose((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_maxrms_and_nbarlayer13_fe_xtrl_loose[0]);
        TGraph gr_acceptance_track_selection_fe_xtrl_loose((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_track_selection_fe_xtrl_loose[0]);
        TGraph gr_acceptance_psd_stk_match_fe_xtrl_loose((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_psd_stk_match_fe_xtrl_loose[0]);
        TGraph gr_acceptance_psd_charge_fe_xtrl_loose((int)gr_acc_bgo_fiducial_energy.size(), &gr_acc_bgo_fiducial_energy[0], &acceptance_psd_charge_fe_xtrl_loose[0]);

        gr_acceptance_bgo_fiducial_het_fe_xtrl_tight.SetName("gr_acceptance_bgo_fiducial_het_fe_xtrl_tight");
        gr_acceptance_nbarlayer13_fe_xtrl_tight.SetName("gr_acceptance_nbarlayer13_fe_xtrl_tight");
        gr_accepatance_maxrms_fe_xtrl_tight.SetName("gr_accepatance_maxrms_fe_xtrl_tight");
        gr_acceptance_maxrms_and_nbarlayer13_fe_xtrl_tight.SetName("gr_acceptance_maxrms_and_nbarlayer13_fe_xtrl_tight");
        gr_acceptance_track_selection_fe_xtrl_tight.SetName("gr_acceptance_track_selection_fe_xtrl_tight");
        gr_acceptance_psd_stk_match_fe_xtrl_tight.SetName("gr_acceptance_psd_stk_match_fe_xtrl_tight");
        gr_acceptance_psd_charge_fe_xtrl_tight.SetName("gr_acceptance_psd_charge_fe_xtrl_tight");

        gr_acceptance_bgo_fiducial_het_fe_xtrl_loose.SetName("gr_acceptance_bgo_fiducial_het_fe_xtrl_loose");
        gr_acceptance_nbarlayer13_fe_xtrl_loose.SetName("gr_acceptance_nbarlayer13_fe_xtrl_loose");
        gr_accepatance_maxrms_fe_xtrl_loose.SetName("gr_accepatance_maxrms_fe_xtrl_loose");
        gr_acceptance_maxrms_and_nbarlayer13_fe_xtrl_loose.SetName("gr_acceptance_maxrms_and_nbarlayer13_fe_xtrl_loose");
        gr_acceptance_track_selection_fe_xtrl_loose.SetName("gr_acceptance_track_selection_fe_xtrl_loose");
        gr_acceptance_psd_stk_match_fe_xtrl_loose.SetName("gr_acceptance_psd_stk_match_fe_xtrl_loose");
        gr_acceptance_psd_charge_fe_xtrl_loose.SetName("gr_acceptance_psd_charge_fe_xtrl_loose");

        // Write output file
        TFile *outfile = TFile::Open(output_file, "RECREATE");
        if (!outfile->IsOpen()) {
            std::cerr << "\n\nError writing output ROOT file [" << output_file << "]\n\n";
            exit(100);
        }

        gr_acceptance_bgo_fiducial_het_fe_xtrl_tight.Write();
        gr_acceptance_nbarlayer13_fe_xtrl_tight.Write();
        gr_accepatance_maxrms_fe_xtrl_tight.Write();
        gr_acceptance_maxrms_and_nbarlayer13_fe_xtrl_tight.Write();
        gr_acceptance_track_selection_fe_xtrl_tight.Write();
        gr_acceptance_psd_stk_match_fe_xtrl_tight.Write();
        gr_acceptance_psd_charge_fe_xtrl_tight.Write();

        gr_acceptance_bgo_fiducial_het_fe_xtrl_loose.Write();
        gr_acceptance_nbarlayer13_fe_xtrl_loose.Write();
        gr_accepatance_maxrms_fe_xtrl_loose.Write();
        gr_acceptance_maxrms_and_nbarlayer13_fe_xtrl_loose.Write();
        gr_acceptance_track_selection_fe_xtrl_loose.Write();
        gr_acceptance_psd_stk_match_fe_xtrl_loose.Write();
        gr_acceptance_psd_charge_fe_xtrl_loose.Write();

        outfile->Close();

    }