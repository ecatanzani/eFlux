#include <iostream>

#include "TH1D.h"
#include "TFile.h"
#include "TEfficiency.h"

void buildMCEfficiency(const char* input_file_path, const char* output_file_path);
void buildDATAEfficiency(const char* input_file_path, const char* output_file_path);

void buildMCEfficiency(const char* input_file_path, const char* output_file_path)
{
    TFile *input_file = TFile::Open(input_file_path, "READ");
    if (!input_file->IsOpen())
    {
        std::cerr << "\n\nError opening input ROOT file [" << input_file_path << "]\n\n";
        exit(100);
    }
    
    auto h_geometric = static_cast<TH1D*>(input_file->Get("h_geometric"));
    auto h_geometric_trigger = static_cast<TH1D*>(input_file->Get("h_geometric_trigger"));
    auto h_trigger = static_cast<TH1D*>(input_file->Get("h_trigger"));
    auto h_maxElayer = static_cast<TH1D*>(input_file->Get("h_maxElayer"));
    auto h_maxBarlayer = static_cast<TH1D*>(input_file->Get("h_maxBarlayer"));
    auto h_BGOTrackContainment = static_cast<TH1D*>(input_file->Get("h_BGOTrackContainment"));
    auto h_bgo_fiducial = static_cast<TH1D*>(input_file->Get("h_bgo_fiducial"));
    auto h_nbarlayer13 = static_cast<TH1D*>(input_file->Get("h_nbarlayer13"));
    auto h_maxrms = static_cast<TH1D*>(input_file->Get("h_maxrms"));
    auto h_track_selection = static_cast<TH1D*>(input_file->Get("h_track_selection"));
    auto h_psd_stk_match = static_cast<TH1D*>(input_file->Get("h_psd_stk_match"));
    auto h_psd_charge = static_cast<TH1D*>(input_file->Get("h_psd_charge"));
    auto h_stk_charge = static_cast<TH1D*>(input_file->Get("h_stk_charge"));
    auto h_all_cuts = static_cast<TH1D*>(input_file->Get("h_all_cuts"));

    h_geometric->SetDirectory(0);
    h_geometric_trigger->SetDirectory(0);
    h_trigger->SetDirectory(0);
    h_maxElayer->SetDirectory(0);
    h_maxBarlayer->SetDirectory(0);
    h_BGOTrackContainment->SetDirectory(0);
    h_bgo_fiducial->SetDirectory(0);
    h_nbarlayer13->SetDirectory(0);
    h_maxrms->SetDirectory(0);
    h_track_selection->SetDirectory(0);
    h_psd_stk_match->SetDirectory(0);
    h_psd_charge->SetDirectory(0);
    h_stk_charge->SetDirectory(0);
    h_all_cuts->SetDirectory(0);

    input_file->Close();

    std::unique_ptr<TEfficiency> trigger_eff;
    std::unique_ptr<TEfficiency> maxElayer_eff;
    std::unique_ptr<TEfficiency> maxBarLayer_eff;
    std::unique_ptr<TEfficiency> BGOTrackContainment_eff;
    std::unique_ptr<TEfficiency> bgo_fiducial_eff;
    std::unique_ptr<TEfficiency> nbarlayer13_eff;
    std::unique_ptr<TEfficiency> maxrms_eff;
    std::unique_ptr<TEfficiency> track_selection_eff;
    std::unique_ptr<TEfficiency> psd_stk_match_eff;
    std::unique_ptr<TEfficiency> psd_charge_eff;
    std::unique_ptr<TEfficiency> stk_charge_eff;
    std::unique_ptr<TEfficiency> all_cuts_eff;

    if (TEfficiency::CheckConsistency(*h_geometric_trigger, *h_geometric))
        trigger_eff = std::make_unique<TEfficiency>(*h_geometric_trigger, *h_geometric);
    if (TEfficiency::CheckConsistency(*h_maxElayer, *h_trigger))
        maxElayer_eff = std::make_unique<TEfficiency>(*h_maxElayer, *h_trigger);
    if (TEfficiency::CheckConsistency(*h_maxBarlayer, *h_trigger))
        maxBarLayer_eff = std::make_unique<TEfficiency>(*h_maxBarlayer, *h_trigger);
    if (TEfficiency::CheckConsistency(*h_BGOTrackContainment, *h_trigger))
        BGOTrackContainment_eff = std::make_unique<TEfficiency>(*h_BGOTrackContainment, *h_trigger);
    if (TEfficiency::CheckConsistency(*h_bgo_fiducial, *h_trigger))
        bgo_fiducial_eff = std::make_unique<TEfficiency>(*h_bgo_fiducial, *h_trigger);
    if (TEfficiency::CheckConsistency(*h_nbarlayer13, *h_bgo_fiducial))
        nbarlayer13_eff = std::make_unique<TEfficiency>(*h_nbarlayer13, *h_bgo_fiducial);
    if (TEfficiency::CheckConsistency(*h_maxrms, *h_nbarlayer13))
        maxrms_eff = std::make_unique<TEfficiency>(*h_maxrms, *h_nbarlayer13);
    if (TEfficiency::CheckConsistency(*h_track_selection, *h_maxrms))
        track_selection_eff = std::make_unique<TEfficiency>(*h_track_selection, *h_maxrms);
    if (TEfficiency::CheckConsistency(*h_psd_stk_match, *h_track_selection))
        psd_stk_match_eff = std::make_unique<TEfficiency>(*h_psd_stk_match, *h_track_selection);
    if (TEfficiency::CheckConsistency(*h_psd_charge, *h_psd_stk_match))
        psd_charge_eff = std::make_unique<TEfficiency>(*h_psd_charge, *h_psd_stk_match);
    if (TEfficiency::CheckConsistency(*h_stk_charge, *h_psd_stk_match))
        stk_charge_eff = std::make_unique<TEfficiency>(*h_stk_charge, *h_psd_stk_match);
    if (TEfficiency::CheckConsistency(*h_all_cuts, *h_trigger))
        all_cuts_eff = std::make_unique<TEfficiency>(*h_all_cuts, *h_trigger);

    trigger_eff->SetStatisticOption(TEfficiency::kBUniform);
    maxElayer_eff->SetStatisticOption(TEfficiency::kBUniform);
    maxBarLayer_eff->SetStatisticOption(TEfficiency::kBUniform);
    BGOTrackContainment_eff->SetStatisticOption(TEfficiency::kBUniform);
    bgo_fiducial_eff->SetStatisticOption(TEfficiency::kBUniform);
    nbarlayer13_eff->SetStatisticOption(TEfficiency::kBUniform);
    maxrms_eff->SetStatisticOption(TEfficiency::kBUniform);
    track_selection_eff->SetStatisticOption(TEfficiency::kBUniform);
    psd_stk_match_eff->SetStatisticOption(TEfficiency::kBUniform);
    psd_charge_eff->SetStatisticOption(TEfficiency::kBUniform);
    stk_charge_eff->SetStatisticOption(TEfficiency::kBUniform);
    all_cuts_eff->SetStatisticOption(TEfficiency::kBUniform);
    
    trigger_eff->SetName("trigger_eff");
    maxElayer_eff->SetName("maxElayer_eff");
    maxBarLayer_eff->SetName("maxBarLayer_eff");
    BGOTrackContainment_eff->SetName("BGOTrackContainment_eff");
    bgo_fiducial_eff->SetName("bgo_fiducial_eff");
    nbarlayer13_eff->SetName("nbarlayer13_eff");
    maxrms_eff->SetName("maxrms_eff");
    track_selection_eff->SetName("track_selection_eff");
    psd_stk_match_eff->SetName("psd_stk_match_eff");
    psd_charge_eff->SetName("psd_charge_eff");
    stk_charge_eff->SetName("stk_charge_eff");
    all_cuts_eff->SetName("all_cuts_eff");

    trigger_eff->SetTitle("Trigger efficiency");
    maxElayer_eff->SetTitle("maxElayer efficiency");
    maxBarLayer_eff->SetTitle("maxBarLayer efficiency");
    BGOTrackContainment_eff->SetTitle("BGO track containment efficiency");
    bgo_fiducial_eff->SetTitle("BGO fiducial volume efficiency");
    nbarlayer13_eff->SetTitle("nbarlayer13 efficiency");
    maxrms_eff->SetTitle("maxrms efficiency");
    track_selection_eff->SetTitle("Track selection efficiency");
    psd_stk_match_eff->SetTitle("PSD/STK match efficiency");
    psd_charge_eff->SetTitle("PSD charge efficiency");
    stk_charge_eff->SetTitle("STK charge efficiency");
    all_cuts_eff->SetTitle("All cuts efficiency");
    
    TFile *output_file = TFile::Open(output_file_path, "RECREATE");
    if (!output_file->IsOpen())
    {
        std::cerr << "\n\nError opening output ROOT file [" << output_file_path << "]\n\n";
        exit(100);
    }

    trigger_eff->Write();
    maxElayer_eff->Write();
    maxBarLayer_eff->Write();
    BGOTrackContainment_eff->Write();
    bgo_fiducial_eff->Write();
    nbarlayer13_eff->Write();
    maxrms_eff->Write();
    track_selection_eff->Write();
    psd_stk_match_eff->Write();
    psd_charge_eff->Write();
    stk_charge_eff->Write();
    all_cuts_eff->Write();
    
    output_file->Write();
    output_file->Close();
}

void buildDATAEfficiency(const char* input_file_path, const char* output_file_path)
{
    TFile *input_file = TFile::Open(input_file_path, "READ");
    if (!input_file->IsOpen())
    {
        std::cerr << "\n\nError opening input ROOT file [" << input_file_path << "]\n\n";
        exit(100);
    }
    
    auto h_trigger = static_cast<TH1D*>(input_file->Get("h_trigger"));
    auto h_maxElayer = static_cast<TH1D*>(input_file->Get("h_maxElayer"));
    auto h_maxBarlayer = static_cast<TH1D*>(input_file->Get("h_maxBarlayer"));
    auto h_BGOTrackContainment = static_cast<TH1D*>(input_file->Get("h_BGOTrackContainment"));
    auto h_bgo_fiducial = static_cast<TH1D*>(input_file->Get("h_bgo_fiducial"));
    auto h_nbarlayer13 = static_cast<TH1D*>(input_file->Get("h_nbarlayer13"));
    auto h_maxrms = static_cast<TH1D*>(input_file->Get("h_maxrms"));
    auto h_track_selection = static_cast<TH1D*>(input_file->Get("h_track_selection"));
    auto h_psd_stk_match = static_cast<TH1D*>(input_file->Get("h_psd_stk_match"));
    auto h_psd_charge = static_cast<TH1D*>(input_file->Get("h_psd_charge"));
    auto h_stk_charge = static_cast<TH1D*>(input_file->Get("h_stk_charge"));
    auto h_all_cuts = static_cast<TH1D*>(input_file->Get("h_all_cuts"));
    
    h_trigger->SetDirectory(0);
    h_maxElayer->SetDirectory(0);
    h_maxBarlayer->SetDirectory(0);
    h_BGOTrackContainment->SetDirectory(0);
    h_bgo_fiducial->SetDirectory(0);
    h_nbarlayer13->SetDirectory(0);
    h_maxrms->SetDirectory(0);
    h_track_selection->SetDirectory(0);
    h_psd_stk_match->SetDirectory(0);
    h_psd_charge->SetDirectory(0);
    h_stk_charge->SetDirectory(0);
    h_all_cuts->SetDirectory(0);

    input_file->Close();

    std::unique_ptr<TEfficiency> maxElayer_eff;
    std::unique_ptr<TEfficiency> maxBarLayer_eff;
    std::unique_ptr<TEfficiency> BGOTrackContainment_eff;
    std::unique_ptr<TEfficiency> bgo_fiducial_eff;
    std::unique_ptr<TEfficiency> nbarlayer13_eff;
    std::unique_ptr<TEfficiency> maxrms_eff;
    std::unique_ptr<TEfficiency> track_selection_eff;
    std::unique_ptr<TEfficiency> psd_stk_match_eff;
    std::unique_ptr<TEfficiency> psd_charge_eff;
    std::unique_ptr<TEfficiency> stk_charge_eff;
    std::unique_ptr<TEfficiency> all_cuts_eff;
    
    if (TEfficiency::CheckConsistency(*h_maxElayer, *h_trigger))
        maxElayer_eff = std::make_unique<TEfficiency>(*h_maxElayer, *h_trigger);
    if (TEfficiency::CheckConsistency(*h_maxBarlayer, *h_trigger))
        maxBarLayer_eff = std::make_unique<TEfficiency>(*h_maxBarlayer, *h_trigger);
    if (TEfficiency::CheckConsistency(*h_BGOTrackContainment, *h_trigger))
        BGOTrackContainment_eff = std::make_unique<TEfficiency>(*h_BGOTrackContainment, *h_trigger);
    if (TEfficiency::CheckConsistency(*h_bgo_fiducial, *h_trigger))
        bgo_fiducial_eff = std::make_unique<TEfficiency>(*h_bgo_fiducial, *h_trigger);
    if (TEfficiency::CheckConsistency(*h_nbarlayer13, *h_bgo_fiducial))
        nbarlayer13_eff = std::make_unique<TEfficiency>(*h_nbarlayer13, *h_bgo_fiducial);
    if (TEfficiency::CheckConsistency(*h_maxrms, *h_nbarlayer13))
        maxrms_eff = std::make_unique<TEfficiency>(*h_maxrms, *h_nbarlayer13);
    if (TEfficiency::CheckConsistency(*h_track_selection, *h_maxrms))
        track_selection_eff = std::make_unique<TEfficiency>(*h_track_selection, *h_maxrms);
    if (TEfficiency::CheckConsistency(*h_psd_stk_match, *h_track_selection))
        psd_stk_match_eff = std::make_unique<TEfficiency>(*h_psd_stk_match, *h_track_selection);
    if (TEfficiency::CheckConsistency(*h_psd_charge, *h_psd_stk_match))
        psd_charge_eff = std::make_unique<TEfficiency>(*h_psd_charge, *h_psd_stk_match);
    if (TEfficiency::CheckConsistency(*h_stk_charge, *h_psd_stk_match))
        stk_charge_eff = std::make_unique<TEfficiency>(*h_stk_charge, *h_psd_stk_match);
    if (TEfficiency::CheckConsistency(*h_all_cuts, *h_trigger))
        all_cuts_eff = std::make_unique<TEfficiency>(*h_all_cuts, *h_trigger);

    maxElayer_eff->SetStatisticOption(TEfficiency::kBUniform);
    maxBarLayer_eff->SetStatisticOption(TEfficiency::kBUniform);
    BGOTrackContainment_eff->SetStatisticOption(TEfficiency::kBUniform);
    bgo_fiducial_eff->SetStatisticOption(TEfficiency::kBUniform);
    nbarlayer13_eff->SetStatisticOption(TEfficiency::kBUniform);
    maxrms_eff->SetStatisticOption(TEfficiency::kBUniform);
    track_selection_eff->SetStatisticOption(TEfficiency::kBUniform);
    psd_stk_match_eff->SetStatisticOption(TEfficiency::kBUniform);
    psd_charge_eff->SetStatisticOption(TEfficiency::kBUniform);
    stk_charge_eff->SetStatisticOption(TEfficiency::kBUniform);
    all_cuts_eff->SetStatisticOption(TEfficiency::kBUniform);

    maxElayer_eff->SetName("maxElayer_eff");
    maxBarLayer_eff->SetName("maxBarLayer_eff");
    BGOTrackContainment_eff->SetName("BGOTrackContainment_eff");
    bgo_fiducial_eff->SetName("bgo_fiducial_eff");
    nbarlayer13_eff->SetName("nbarlayer13_eff");
    maxrms_eff->SetName("maxrms_eff");
    track_selection_eff->SetName("track_selection_eff");
    psd_stk_match_eff->SetName("psd_stk_match_eff");
    psd_charge_eff->SetName("psd_charge_eff");
    stk_charge_eff->SetName("stk_charge_eff");
    all_cuts_eff->SetName("all_cuts_eff");
    
    maxElayer_eff->SetTitle("maxElayer efficiency");
    maxBarLayer_eff->SetTitle("maxBarLayer efficiency");
    BGOTrackContainment_eff->SetTitle("BGO track containment efficiency");
    bgo_fiducial_eff->SetTitle("BGO fiducial volume efficiency");
    nbarlayer13_eff->SetTitle("nbarlayer13 efficiency");
    maxrms_eff->SetTitle("maxrms efficiency");
    track_selection_eff->SetTitle("Track selection efficiency");
    psd_stk_match_eff->SetTitle("PSD/STK match efficiency");
    psd_charge_eff->SetTitle("PSD charge efficiency");
    stk_charge_eff->SetTitle("STK charge efficiency");
    all_cuts_eff->SetTitle("All cuts efficiency");
    
    TFile *output_file = TFile::Open(output_file_path, "RECREATE");
    if (!output_file->IsOpen())
    {
        std::cerr << "\n\nError opening output ROOT file [" << output_file_path << "]\n\n";
        exit(100);
    }
    
    maxElayer_eff->Write();
    maxBarLayer_eff->Write();
    BGOTrackContainment_eff->Write();
    bgo_fiducial_eff->Write();
    nbarlayer13_eff->Write();
    maxrms_eff->Write();
    track_selection_eff->Write();
    psd_stk_match_eff->Write();
    psd_charge_eff->Write();
    stk_charge_eff->Write();
    all_cuts_eff->Write();
    
    output_file->Write();
    output_file->Close();
}