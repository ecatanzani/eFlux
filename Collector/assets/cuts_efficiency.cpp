#include <memory>
#include <iostream>

#include "TH1D.h"
#include "TFile.h"
#include "TEfficiency.h"

void buildEfficiency(const char* input_file_path, const char* output_file_path) {
    
    TFile *input_file = TFile::Open(input_file_path, "READ");
    if (!input_file->IsOpen()) {
        std::cerr << "\n\nError opening input ROOT file [" << input_file_path << "]\n\n";
        exit(100);
    }

    auto h_maxelayer_lastcut_pass = static_cast<TH1D*>(input_file->Get("Stats/h_maxelayer_lastcut_pass"));
    auto h_maxelayer_lastcut = static_cast<TH1D*>(input_file->Get("Stats/h_maxelayer_lastcut"));
    auto h_maxbarlayer_lastcut_pass = static_cast<TH1D*>(input_file->Get("Stats/h_maxbarlayer_lastcut_pass"));
    auto h_maxbarlayer_lastcut = static_cast<TH1D*>(input_file->Get("Stats/h_maxbarlayer_lastcut"));
    auto h_bgotrack_lastcut_pass = static_cast<TH1D*>(input_file->Get("Stats/h_bgotrack_lastcut_pass"));
    auto h_bgotrack_lastcut = static_cast<TH1D*>(input_file->Get("Stats/h_bgotrack_lastcut"));
    auto h_bgofiducial_lastcut_pass = static_cast<TH1D*>(input_file->Get("Stats/h_bgofiducial_lastcut_pass"));
    auto h_bgofiducial_lastcut = static_cast<TH1D*>(input_file->Get("Stats/h_bgofiducial_lastcut"));
    auto h_nbarlayer13_lastcut_pass = static_cast<TH1D*>(input_file->Get("Stats/h_nbarlayer13_lastcut_pass"));
    auto h_nbarlayer13_lastcut = static_cast<TH1D*>(input_file->Get("Stats/h_nbarlayer13_lastcut"));
    auto h_maxrms_lastcut_pass = static_cast<TH1D*>(input_file->Get("Stats/h_maxrms_lastcut_pass"));
    auto h_maxrms_lastcut = static_cast<TH1D*>(input_file->Get("Stats/h_maxrms_lastcut"));
    auto h_trackselection_lastcut_pass = static_cast<TH1D*>(input_file->Get("Stats/h_trackselection_lastcut_pass"));
    auto h_trackselection_lastcut = static_cast<TH1D*>(input_file->Get("Stats/h_trackselection_lastcut"));
    auto h_psdstkmatch_lastcut_pass = static_cast<TH1D*>(input_file->Get("Stats/h_psdstkmatch_lastcut_pass"));
    auto h_psdstkmatch_lastcut = static_cast<TH1D*>(input_file->Get("Stats/h_psdstkmatch_lastcut"));
    auto h_psdcharge_lastcut_pass = static_cast<TH1D*>(input_file->Get("Stats/h_psdcharge_lastcut_pass"));
    auto h_psdcharge_lastcut = static_cast<TH1D*>(input_file->Get("Stats/h_psdcharge_lastcut"));
    auto h_stkcharge_lastcut_pass = static_cast<TH1D*>(input_file->Get("Stats/h_stkcharge_lastcut_pass"));
    auto h_stkcharge_lastcut = static_cast<TH1D*>(input_file->Get("Stats/h_stkcharge_lastcut"));

    h_maxelayer_lastcut_pass->SetDirectory(0);
    h_maxelayer_lastcut->SetDirectory(0);
    h_maxbarlayer_lastcut_pass->SetDirectory(0);
    h_maxbarlayer_lastcut->SetDirectory(0);
    h_bgotrack_lastcut_pass->SetDirectory(0);
    h_bgotrack_lastcut->SetDirectory(0);
    h_bgofiducial_lastcut_pass->SetDirectory(0);
    h_bgofiducial_lastcut->SetDirectory(0);
    h_nbarlayer13_lastcut_pass->SetDirectory(0);
    h_nbarlayer13_lastcut->SetDirectory(0);
    h_maxrms_lastcut_pass->SetDirectory(0);
    h_maxrms_lastcut->SetDirectory(0);
    h_trackselection_lastcut_pass->SetDirectory(0);
    h_trackselection_lastcut->SetDirectory(0);
    h_psdstkmatch_lastcut_pass->SetDirectory(0);
    h_psdstkmatch_lastcut->SetDirectory(0);
    h_psdcharge_lastcut_pass->SetDirectory(0);
    h_psdcharge_lastcut->SetDirectory(0);
    h_stkcharge_lastcut_pass->SetDirectory(0);
    h_stkcharge_lastcut->SetDirectory(0);

    input_file->Close();

    std::unique_ptr<TEfficiency> maxelayer_eff;
    std::unique_ptr<TEfficiency> maxbarlayer_eff;
    std::unique_ptr<TEfficiency> bgotrack_eff;
    std::unique_ptr<TEfficiency> bgofiducial_eff;
    std::unique_ptr<TEfficiency> nbarlayer13_eff;
    std::unique_ptr<TEfficiency> maxrms_eff;
    std::unique_ptr<TEfficiency> trackselection_eff;
    std::unique_ptr<TEfficiency> psdstkmatch_eff;
    std::unique_ptr<TEfficiency> psdcharge_eff;
    std::unique_ptr<TEfficiency> stkcharge_eff;

    if (TEfficiency::CheckConsistency(*h_maxelayer_lastcut_pass, *h_maxelayer_lastcut))
        maxelayer_eff = std::make_unique<TEfficiency>(*h_maxelayer_lastcut_pass, *h_maxelayer_lastcut);
    if (TEfficiency::CheckConsistency(*h_maxbarlayer_lastcut_pass, *h_maxbarlayer_lastcut))
        maxbarlayer_eff = std::make_unique<TEfficiency>(*h_maxbarlayer_lastcut_pass, *h_maxbarlayer_lastcut);
    if (TEfficiency::CheckConsistency(*h_bgotrack_lastcut_pass, *h_bgotrack_lastcut))
        bgotrack_eff = std::make_unique<TEfficiency>(*h_bgotrack_lastcut_pass, *h_bgotrack_lastcut);
    if (TEfficiency::CheckConsistency(*h_bgofiducial_lastcut_pass, *h_bgofiducial_lastcut))
        bgofiducial_eff = std::make_unique<TEfficiency>(*h_bgofiducial_lastcut_pass, *h_bgofiducial_lastcut);
    if (TEfficiency::CheckConsistency(*h_nbarlayer13_lastcut_pass, *h_nbarlayer13_lastcut))
        nbarlayer13_eff = std::make_unique<TEfficiency>(*h_nbarlayer13_lastcut_pass, *h_nbarlayer13_lastcut);
    if (TEfficiency::CheckConsistency(*h_maxrms_lastcut_pass, *h_maxrms_lastcut))
        maxrms_eff = std::make_unique<TEfficiency>(*h_maxrms_lastcut_pass, *h_maxrms_lastcut);
    if (TEfficiency::CheckConsistency(*h_trackselection_lastcut_pass, *h_trackselection_lastcut))
        trackselection_eff = std::make_unique<TEfficiency>(*h_trackselection_lastcut_pass, *h_trackselection_lastcut);
    if (TEfficiency::CheckConsistency(*h_psdstkmatch_lastcut_pass, *h_psdstkmatch_lastcut))
        psdstkmatch_eff = std::make_unique<TEfficiency>(*h_psdstkmatch_lastcut_pass, *h_psdstkmatch_lastcut);
    if (TEfficiency::CheckConsistency(*h_psdcharge_lastcut_pass, *h_psdcharge_lastcut))
        psdcharge_eff = std::make_unique<TEfficiency>(*h_psdcharge_lastcut_pass, *h_psdcharge_lastcut);
    if (TEfficiency::CheckConsistency(*h_stkcharge_lastcut_pass, *h_stkcharge_lastcut))
        stkcharge_eff = std::make_unique<TEfficiency>(*h_stkcharge_lastcut_pass, *h_stkcharge_lastcut);

    maxelayer_eff->SetStatisticOption(TEfficiency::kBUniform);
    maxbarlayer_eff->SetStatisticOption(TEfficiency::kBUniform);
    bgotrack_eff->SetStatisticOption(TEfficiency::kBUniform);
    bgofiducial_eff->SetStatisticOption(TEfficiency::kBUniform);
    nbarlayer13_eff->SetStatisticOption(TEfficiency::kBUniform);
    maxrms_eff->SetStatisticOption(TEfficiency::kBUniform);
    trackselection_eff->SetStatisticOption(TEfficiency::kBUniform);
    psdstkmatch_eff->SetStatisticOption(TEfficiency::kBUniform);
    psdcharge_eff->SetStatisticOption(TEfficiency::kBUniform);
    stkcharge_eff->SetStatisticOption(TEfficiency::kBUniform);

    maxelayer_eff->SetName("maxelayer_eff");
    maxbarlayer_eff->SetName("maxbarlayer_eff");
    bgotrack_eff->SetName("bgotrack_eff");
    bgofiducial_eff->SetName("bgofiducial_eff");
    nbarlayer13_eff->SetName("nbarlayer13_eff");
    maxrms_eff->SetName("maxrms_eff");
    trackselection_eff->SetName("trackselection_eff");
    psdstkmatch_eff->SetName("psdstkmatch_eff");
    psdcharge_eff->SetName("psdcharge_eff");
    stkcharge_eff->SetName("stkcharge_eff");

    maxelayer_eff->SetTitle("MaxELayer efficiency");
    maxbarlayer_eff->SetTitle("MaxBarLayer efficiency");
    bgotrack_eff->SetTitle("BGO Track efficiency");
    bgofiducial_eff->SetTitle("BGO Fiducial efficiency");
    nbarlayer13_eff->SetTitle("nBarLayer13 efficiency");
    maxrms_eff->SetTitle("max RMS efficiency");
    trackselection_eff->SetTitle("Track Selection efficiency");
    psdstkmatch_eff->SetTitle("PSD/STK Match efficiency");
    psdcharge_eff->SetTitle("PSD charge efficiency");
    stkcharge_eff->SetTitle("STK charge efficiency");

    TFile *output_file = TFile::Open(output_file_path, "RECREATE");
    if (!output_file->IsOpen()) {
        std::cerr << "\n\nError opening output ROOT file [" << output_file_path << "]\n\n";
        exit(100);
    }

    maxelayer_eff->Write();
    maxbarlayer_eff->Write();
    bgotrack_eff->Write();
    bgofiducial_eff->Write();
    nbarlayer13_eff->Write();
    maxrms_eff->Write();
    trackselection_eff->Write();
    psdstkmatch_eff->Write();
    psdcharge_eff->Write();
    stkcharge_eff->Write();

    output_file->Close();
}