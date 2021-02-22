#include <iostream>

#include "TH1D.h"
#include "TFile.h"
#include "TEfficiency.h"

void BuildEfficiency(
    const char* input_file_path, 
    const char* output_file_path="efficiencyout.root",
    const bool mc=false)
{
    TFile *input_file = TFile::Open(input_file_path, "READ");
    if (!input_file->IsOpen())
    {
        std::cerr << "\n\nError opening input ROOT file [" << input_file_path << "]\n\n";
        exit(100);
    }

    // Read histos
    /*
    TH1D* h_geo_before_trigger_cut = nullptr;
    if (mc)
        h_geo_before_trigger_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_geo_before_trigger_cut"));
    */
    auto h_trigger_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_trigger_cut"));
    auto h_geometric_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_geometric_cut"));
    auto h_maxElayer_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_maxElayer_cut"));
    auto h_maxBarlayer_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_maxBarlayer_cut"));
    auto h_BGOTrackContainment_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_BGOTrackContainment_cut"));
    auto h_bgo_fiducial_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_bgo_fiducial_cut"));
    auto h_nbarlayer13_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_nbarlayer13_cut"));
    auto h_maxrms_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_maxrms_cut"));
    auto h_track_selection_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_track_selection_cut"));
    auto h_psd_stk_match_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_psd_stk_match_cut"));
    auto h_psd_charge_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_psd_charge_cut"));
    auto h_stk_charge_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_stk_charge_cut"));
    auto h_all_cuts_cut = static_cast<TH1D*>(input_file->Get("Cuts/h_all_cuts_cut"));
    
    // Unlink histos
    /*
    if (mc)
        h_geo_before_trigger_cut->SetDirectory(0);
    */
    h_trigger_cut->SetDirectory(0);
    h_geometric_cut->SetDirectory(0);
    h_maxElayer_cut->SetDirectory(0);
    h_maxBarlayer_cut->SetDirectory(0);
    h_BGOTrackContainment_cut->SetDirectory(0);
    h_bgo_fiducial_cut->SetDirectory(0);
    h_nbarlayer13_cut->SetDirectory(0);
    h_maxrms_cut->SetDirectory(0);
    h_track_selection_cut->SetDirectory(0);
    h_psd_stk_match_cut->SetDirectory(0);
    h_psd_charge_cut->SetDirectory(0);
    h_stk_charge_cut->SetDirectory(0);
    h_all_cuts_cut->SetDirectory(0);

    input_file->Close();

    // Build efficiencies
    //std::unique_ptr<TEfficiency> trigger_eff;
    std::unique_ptr<TEfficiency> geometric_factor_eff;
    std::unique_ptr<TEfficiency> maxElayer_cut_eff;
    std::unique_ptr<TEfficiency> maxBarLayer_cut_eff;
    std::unique_ptr<TEfficiency> BGOTrackContainment_cut_eff;
    std::unique_ptr<TEfficiency> bgo_fiducial_cut_eff;
    std::unique_ptr<TEfficiency> nbarlayer13_cut_eff;
    std::unique_ptr<TEfficiency> maxrms_cut_eff;
    std::unique_ptr<TEfficiency> track_selection_eff;
    std::unique_ptr<TEfficiency> psd_stk_match_eff;
    std::unique_ptr<TEfficiency> psd_charge_eff;
    std::unique_ptr<TEfficiency> stk_charge_eff;
    std::unique_ptr<TEfficiency> all_cuts_eff;

    /*
    if (mc)
    {
        if (TEfficiency::CheckConsistency(*h_trigger_cut, *h_geo_before_trigger_cut))
            trigger_eff = std::make_unique<TEfficiency>(*h_trigger_cut, *h_geo_before_trigger_cut);
    }
    */
    if (TEfficiency::CheckConsistency(*h_geometric_cut, *h_trigger_cut))
        geometric_factor_eff = std::make_unique<TEfficiency>(*h_geometric_cut, *h_trigger_cut);
    if (TEfficiency::CheckConsistency(*h_maxElayer_cut, *h_trigger_cut))
        maxElayer_cut_eff = std::make_unique<TEfficiency>(*h_maxElayer_cut, *h_trigger_cut);
    if (TEfficiency::CheckConsistency(*h_maxBarlayer_cut, *h_trigger_cut))
        maxBarLayer_cut_eff = std::make_unique<TEfficiency>(*h_maxBarlayer_cut, *h_trigger_cut);
    if (TEfficiency::CheckConsistency(*h_BGOTrackContainment_cut, *h_trigger_cut))
        BGOTrackContainment_cut_eff = std::make_unique<TEfficiency>(*h_BGOTrackContainment_cut, *h_trigger_cut);
    if (TEfficiency::CheckConsistency(*h_bgo_fiducial_cut, *h_trigger_cut))
        bgo_fiducial_cut_eff = std::make_unique<TEfficiency>(*h_bgo_fiducial_cut, *h_trigger_cut);
    if (TEfficiency::CheckConsistency(*h_nbarlayer13_cut, *h_bgo_fiducial_cut))
        nbarlayer13_cut_eff = std::make_unique<TEfficiency>(*h_nbarlayer13_cut, *h_bgo_fiducial_cut);
    if (TEfficiency::CheckConsistency(*h_maxrms_cut, *h_nbarlayer13_cut))
        maxrms_cut_eff = std::make_unique<TEfficiency>(*h_maxrms_cut, *h_nbarlayer13_cut);
    if (TEfficiency::CheckConsistency(*h_track_selection_cut, *h_maxrms_cut))
        track_selection_eff = std::make_unique<TEfficiency>(*h_track_selection_cut, *h_maxrms_cut);
    if (TEfficiency::CheckConsistency(*h_psd_stk_match_cut, *h_track_selection_cut))
        psd_stk_match_eff = std::make_unique<TEfficiency>(*h_psd_stk_match_cut, *h_track_selection_cut);
    if (TEfficiency::CheckConsistency(*h_psd_charge_cut, *h_psd_stk_match_cut))
        psd_charge_eff = std::make_unique<TEfficiency>(*h_psd_charge_cut, *h_psd_stk_match_cut);
    if (TEfficiency::CheckConsistency(*h_stk_charge_cut, *h_psd_charge_cut))
        stk_charge_eff = std::make_unique<TEfficiency>(*h_stk_charge_cut, *h_psd_charge_cut);
    if (TEfficiency::CheckConsistency(*h_all_cuts_cut, *h_trigger_cut))
        all_cuts_eff = std::make_unique<TEfficiency>(*h_all_cuts_cut, *h_trigger_cut);

    /*
    if (mc)
        trigger_eff->SetStatisticOption(TEfficiency::kBUniform);
    */
    geometric_factor_eff->SetStatisticOption(TEfficiency::kBUniform);
    maxElayer_cut_eff->SetStatisticOption(TEfficiency::kBUniform);
    maxBarLayer_cut_eff->SetStatisticOption(TEfficiency::kBUniform);
    BGOTrackContainment_cut_eff->SetStatisticOption(TEfficiency::kBUniform);
    bgo_fiducial_cut_eff->SetStatisticOption(TEfficiency::kBUniform);
    nbarlayer13_cut_eff->SetStatisticOption(TEfficiency::kBUniform);
    maxrms_cut_eff->SetStatisticOption(TEfficiency::kBUniform);
    track_selection_eff->SetStatisticOption(TEfficiency::kBUniform);
    psd_stk_match_eff->SetStatisticOption(TEfficiency::kBUniform);
    psd_charge_eff->SetStatisticOption(TEfficiency::kBUniform);
    stk_charge_eff->SetStatisticOption(TEfficiency::kBUniform);
    all_cuts_eff->SetStatisticOption(TEfficiency::kBUniform);

    /*
    if (mc)
        trigger_eff->SetName("trigger_eff");
    */
    geometric_factor_eff->SetName("geometric_factor_eff");
    maxElayer_cut_eff->SetName("maxElayer_cut_eff");
    maxBarLayer_cut_eff->SetName("maxBarLayer_cut_eff");
    BGOTrackContainment_cut_eff->SetName("BGOTrackContainment_cut_eff");
    bgo_fiducial_cut_eff->SetName("bgo_fiducial_cut_eff");
    nbarlayer13_cut_eff->SetName("nbarlayer13_cut_eff");
    maxrms_cut_eff->SetName("maxrms_cut_eff");
    track_selection_eff->SetName("track_selection_eff");
    psd_stk_match_eff->SetName("psd_stk_match_eff");
    psd_charge_eff->SetName("psd_charge_eff");
    stk_charge_eff->SetName("stk_charge_eff");
    all_cuts_eff->SetName("all_cuts_eff");

    /*
    if (mc)
        trigger_eff->SetTitle("Trigger efficiency");
    */
    geometric_factor_eff->SetTitle("Geometric factor efficiency");
    maxElayer_cut_eff->SetTitle("maxElayer cut efficiency");
    maxBarLayer_cut_eff->SetTitle("maxBarLayer cut efficiency");
    BGOTrackContainment_cut_eff->SetTitle("BGOTrackContainment cut efficiency");
    bgo_fiducial_cut_eff->SetTitle("BGO fiducial cut efficiency");
    nbarlayer13_cut_eff->SetTitle("nbarlayer13 cut efficiency");
    maxrms_cut_eff->SetTitle("maxrms cut efficiency");
    track_selection_eff->SetTitle("track selection efficiency");
    psd_stk_match_eff->SetTitle("PSD-STK match efficiency");
    psd_charge_eff->SetTitle("PSD charge efficiency");
    stk_charge_eff->SetTitle("STK charge efficiency");
    all_cuts_eff->SetTitle("all cuts efficiency");

    // Write effi iency histos on disk
    TFile *output_file = TFile::Open(output_file_path, "RECREATE");
    if (!output_file->IsOpen())
    {
        std::cerr << "\n\nError opening output ROOT file [" << output_file_path << "]\n\n";
        exit(100);
    }

    /*
    if (mc)
        trigger_eff->Write();
    */
    geometric_factor_eff->Write();
    maxElayer_cut_eff->Write();
    maxBarLayer_cut_eff->Write();
    BGOTrackContainment_cut_eff->Write();
    bgo_fiducial_cut_eff->Write();
    nbarlayer13_cut_eff->Write();
    maxrms_cut_eff->Write();
    track_selection_eff->Write();
    psd_stk_match_eff->Write();
    psd_charge_eff->Write();
    stk_charge_eff->Write();
    all_cuts_eff->Write();
    
    output_file->Write();
    output_file->Close();
}