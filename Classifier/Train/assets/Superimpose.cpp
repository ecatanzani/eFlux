#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"

void Superimpose(
    const char* infile_data,
    const char* infile_mc,
    const char* outputfile)
{
    TFile* _file_data = TFile::Open(infile_data, "READ");
    if (_file_data->IsZombie())
    {
        std::cerr << "\n\nError opening input file [" << infile_data << "]" << std::endl;
        exit(100);
    }

    auto h_energy_data = static_cast<TH1D*>(_file_data->Get("h_energy"));
    h_energy_data->SetName("h_energy_data");
    h_energy_data->SetDirectory(0);
    h_energy_data->SetLineWidth(3);
    h_energy_data->SetLineColor(kBlue);
    auto h_energy_corr_data = static_cast<TH1D*>(_file_data->Get("h_energy_corr"));
    h_energy_corr_data->SetName("h_energy_corr_data");
    h_energy_corr_data->SetDirectory(0);
    h_energy_corr_data->SetLineWidth(3);
    h_energy_corr_data->SetLineColor(kBlue);
    auto h_sumRms_data = static_cast<TH1D*>(_file_data->Get("h_sumRms"));
    h_sumRms_data->SetName("h_sumRms_data");
    h_sumRms_data->SetDirectory(0);
    h_sumRms_data->SetLineWidth(3);
    h_sumRms_data->SetLineColor(kBlue);
    auto h_fracLast_data = static_cast<TH1D*>(_file_data->Get("h_fracLast"));
    h_fracLast_data->SetName("h_fracLast_data");
    h_fracLast_data->SetDirectory(0);
    h_fracLast_data->SetLineWidth(3);
    h_fracLast_data->SetLineColor(kBlue);

    auto h_rmsLayer_1_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_1"));
    h_rmsLayer_1_data->SetName("h_rmsLayer_1_data");
    h_rmsLayer_1_data->SetDirectory(0);
    h_rmsLayer_1_data->SetLineWidth(3);
    h_rmsLayer_1_data->SetLineColor(kBlue);
    auto h_rmsLayer_2_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_2"));
    h_rmsLayer_2_data->SetName("h_rmsLayer_2_data");
    h_rmsLayer_2_data->SetDirectory(0);
    h_rmsLayer_2_data->SetLineWidth(3);
    h_rmsLayer_2_data->SetLineColor(kBlue);
    auto h_rmsLayer_3_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_3"));
    h_rmsLayer_3_data->SetName("h_rmsLayer_3_data");
    h_rmsLayer_3_data->SetDirectory(0);
    h_rmsLayer_3_data->SetLineWidth(3);
    h_rmsLayer_3_data->SetLineColor(kBlue);
    auto h_rmsLayer_4_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_4"));
    h_rmsLayer_4_data->SetName("h_rmsLayer_4_data");
    h_rmsLayer_4_data->SetDirectory(0);
    h_rmsLayer_4_data->SetLineWidth(3);
    h_rmsLayer_4_data->SetLineColor(kBlue);
    auto h_rmsLayer_5_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_5"));
    h_rmsLayer_5_data->SetName("h_rmsLayer_5_data");
    h_rmsLayer_5_data->SetDirectory(0);
    h_rmsLayer_5_data->SetLineWidth(3);
    h_rmsLayer_5_data->SetLineColor(kBlue);
    auto h_rmsLayer_6_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_6"));
    h_rmsLayer_6_data->SetName("h_rmsLayer_6_data");
    h_rmsLayer_6_data->SetDirectory(0);
    h_rmsLayer_6_data->SetLineWidth(3);
    h_rmsLayer_6_data->SetLineColor(kBlue);
    auto h_rmsLayer_7_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_7"));
    h_rmsLayer_7_data->SetName("h_rmsLayer_7_data");
    h_rmsLayer_7_data->SetDirectory(0);
    h_rmsLayer_7_data->SetLineWidth(3);
    h_rmsLayer_7_data->SetLineColor(kBlue);
    auto h_rmsLayer_8_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_8"));
    h_rmsLayer_8_data->SetName("h_rmsLayer_8_data");
    h_rmsLayer_8_data->SetDirectory(0);
    h_rmsLayer_8_data->SetLineWidth(3);
    h_rmsLayer_8_data->SetLineColor(kBlue);
    auto h_rmsLayer_9_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_9"));
    h_rmsLayer_9_data->SetName("h_rmsLayer_9_data");
    h_rmsLayer_9_data->SetDirectory(0);
    h_rmsLayer_9_data->SetLineWidth(3);
    h_rmsLayer_9_data->SetLineColor(kBlue);
    auto h_rmsLayer_10_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_10"));
    h_rmsLayer_10_data->SetName("h_rmsLayer_10_data");
    h_rmsLayer_10_data->SetDirectory(0);
    h_rmsLayer_10_data->SetLineWidth(3);
    h_rmsLayer_10_data->SetLineColor(kBlue);
    auto h_rmsLayer_11_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_11"));
    h_rmsLayer_11_data->SetName("h_rmsLayer_11_data");
    h_rmsLayer_11_data->SetDirectory(0);
    h_rmsLayer_11_data->SetLineWidth(3);
    h_rmsLayer_11_data->SetLineColor(kBlue);
    auto h_rmsLayer_12_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_12"));
    h_rmsLayer_12_data->SetName("h_rmsLayer_12_data");
    h_rmsLayer_12_data->SetDirectory(0);
    h_rmsLayer_12_data->SetLineWidth(3);
    h_rmsLayer_12_data->SetLineColor(kBlue);
    auto h_rmsLayer_13_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_13"));
    h_rmsLayer_13_data->SetName("h_rmsLayer_13_data");
    h_rmsLayer_13_data->SetDirectory(0);
    h_rmsLayer_13_data->SetLineWidth(3);
    h_rmsLayer_13_data->SetLineColor(kBlue);
    auto h_rmsLayer_14_data = static_cast<TH1D*>(_file_data->Get("h_rmsLayer_14"));
    h_rmsLayer_14_data->SetName("h_rmsLayer_14_data");
    h_rmsLayer_14_data->SetDirectory(0);
    h_rmsLayer_14_data->SetLineWidth(3);
    h_rmsLayer_14_data->SetLineColor(kBlue);

    auto h_fracLayer_1_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_1"));
    h_fracLayer_1_data->SetName("h_fracLayer_1_data");
    h_fracLayer_1_data->SetDirectory(0);
    h_fracLayer_1_data->SetLineWidth(3);
    h_fracLayer_1_data->SetLineColor(kBlue);
    auto h_fracLayer_2_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_2"));
    h_fracLayer_2_data->SetName("h_fracLayer_2_data");
    h_fracLayer_2_data->SetDirectory(0);
    h_fracLayer_2_data->SetLineWidth(3);
    h_fracLayer_2_data->SetLineColor(kBlue);
    auto h_fracLayer_3_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_3"));
    h_fracLayer_3_data->SetName("h_fracLayer_3_data");
    h_fracLayer_3_data->SetDirectory(0);
    h_fracLayer_3_data->SetLineWidth(3);
    h_fracLayer_3_data->SetLineColor(kBlue);
    auto h_fracLayer_4_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_4"));
    h_fracLayer_4_data->SetName("h_fracLayer_4_data");
    h_fracLayer_4_data->SetDirectory(0);
    h_fracLayer_4_data->SetLineWidth(3);
    h_fracLayer_4_data->SetLineColor(kBlue);
    auto h_fracLayer_5_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_5"));
    h_fracLayer_5_data->SetName("h_fracLayer_5_data");
    h_fracLayer_5_data->SetDirectory(0);
    h_fracLayer_5_data->SetLineWidth(3);
    h_fracLayer_5_data->SetLineColor(kBlue);
    auto h_fracLayer_6_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_6"));
    h_fracLayer_6_data->SetName("h_fracLayer_6_data");
    h_fracLayer_6_data->SetDirectory(0);
    h_fracLayer_6_data->SetLineWidth(3);
    h_fracLayer_6_data->SetLineColor(kBlue);
    auto h_fracLayer_7_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_7"));
    h_fracLayer_7_data->SetName("h_fracLayer_7_data");
    h_fracLayer_7_data->SetDirectory(0);
    h_fracLayer_7_data->SetLineWidth(3);
    h_fracLayer_7_data->SetLineColor(kBlue);
    auto h_fracLayer_8_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_8"));
    h_fracLayer_8_data->SetName("h_fracLayer_8_data");
    h_fracLayer_8_data->SetDirectory(0);
    h_fracLayer_8_data->SetLineWidth(3);
    h_fracLayer_8_data->SetLineColor(kBlue);
    auto h_fracLayer_9_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_9"));
    h_fracLayer_9_data->SetName("h_fracLayer_9_data");
    h_fracLayer_9_data->SetDirectory(0);
    h_fracLayer_9_data->SetLineWidth(3);
    h_fracLayer_9_data->SetLineColor(kBlue);
    auto h_fracLayer_10_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_10"));
    h_fracLayer_10_data->SetName("h_fracLayer_10_data");
    h_fracLayer_10_data->SetDirectory(0);
    h_fracLayer_10_data->SetLineWidth(3);
    h_fracLayer_10_data->SetLineColor(kBlue);
    auto h_fracLayer_11_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_11"));
    h_fracLayer_11_data->SetName("h_fracLayer_11_data");
    h_fracLayer_11_data->SetDirectory(0);
    h_fracLayer_11_data->SetLineWidth(3);
    h_fracLayer_11_data->SetLineColor(kBlue);
    auto h_fracLayer_12_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_12"));
    h_fracLayer_12_data->SetName("h_fracLayer_12_data");
    h_fracLayer_12_data->SetDirectory(0);
    h_fracLayer_12_data->SetLineWidth(3);
    h_fracLayer_12_data->SetLineColor(kBlue);
    auto h_fracLayer_13_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_13"));
    h_fracLayer_13_data->SetName("h_fracLayer_13_data");
    h_fracLayer_13_data->SetDirectory(0);
    h_fracLayer_13_data->SetLineWidth(3);
    h_fracLayer_13_data->SetLineColor(kBlue);
    auto h_fracLayer_14_data = static_cast<TH1D*>(_file_data->Get("h_fracLayer_14"));
    h_fracLayer_14_data->SetName("h_fracLayer_14_data");
    h_fracLayer_14_data->SetDirectory(0);
    h_fracLayer_14_data->SetLineWidth(3);
    h_fracLayer_14_data->SetLineColor(kBlue);

    auto h_lastbgolayer_data = static_cast<TH1D*>(_file_data->Get("h_lastbgolayer"));
    h_lastbgolayer_data->SetName("h_lastbgolayer_data");
    h_lastbgolayer_data->SetDirectory(0);
    h_lastbgolayer_data->SetLineWidth(3);
    h_lastbgolayer_data->SetLineColor(kBlue);

    auto h_nbgoentries_data = static_cast<TH1D*>(_file_data->Get("h_nbgoentries"));
    h_nbgoentries_data->SetName("h_nbgoentries_data");
    h_nbgoentries_data->SetDirectory(0);
    h_nbgoentries_data->SetLineWidth(3);
    h_nbgoentries_data->SetLineColor(kBlue);

    auto h_xtrl_data = static_cast<TH1D*>(_file_data->Get("h_xtrl"));
    h_xtrl_data->SetName("h_xtrl_data");
    h_xtrl_data->SetDirectory(0);
    h_xtrl_data->SetLineWidth(3);
    h_xtrl_data->SetLineColor(kBlue);

    auto h_nudadc_1_data = static_cast<TH1D*>(_file_data->Get("h_nudadc_1"));
    h_nudadc_1_data->SetName("h_nudadc_1_data");
    h_nudadc_1_data->SetDirectory(0);
    h_nudadc_1_data->SetLineWidth(3);
    h_nudadc_1_data->SetLineColor(kBlue);
    auto h_nudadc_2_data = static_cast<TH1D*>(_file_data->Get("h_nudadc_2"));
    h_nudadc_2_data->SetName("h_nudadc_2_data");
    h_nudadc_2_data->SetDirectory(0);
    h_nudadc_2_data->SetLineWidth(3);
    h_nudadc_2_data->SetLineColor(kBlue);
    auto h_nudadc_3_data = static_cast<TH1D*>(_file_data->Get("h_nudadc_3"));
    h_nudadc_3_data->SetName("h_nudadc_3_data");
    h_nudadc_3_data->SetDirectory(0);
    h_nudadc_3_data->SetLineWidth(3);
    h_nudadc_3_data->SetLineColor(kBlue);
    auto h_nudadc_4_data = static_cast<TH1D*>(_file_data->Get("h_nudadc_4"));
    h_nudadc_4_data->SetName("h_nudadc_4_data");
    h_nudadc_4_data->SetDirectory(0);
    h_nudadc_4_data->SetLineWidth(3);
    h_nudadc_4_data->SetLineColor(kBlue);

    auto h_nud_totaladc_data = static_cast<TH1D*>(_file_data->Get("h_nud_totaladc"));
    h_nud_totaladc_data->SetName("h_nud_totaladc_data");
    h_nud_totaladc_data->SetDirectory(0);
    h_nud_totaladc_data->SetLineWidth(3);
    h_nud_totaladc_data->SetLineColor(kBlue);
    auto h_nud_maxadc_data = static_cast<TH1D*>(_file_data->Get("h_nud_maxadc"));
    h_nud_maxadc_data->SetName("h_nud_maxadc_data");
    h_nud_maxadc_data->SetDirectory(0);
    h_nud_maxadc_data->SetLineWidth(3);
    h_nud_maxadc_data->SetLineColor(kBlue);

    _file_data->Close();

    TFile* _file_mc = TFile::Open(infile_mc, "READ");
    if (_file_mc->IsZombie())
    {
        std::cerr << "\n\nError opening input file [" << infile_data << "]" << std::endl;
        exit(100);
    }

    auto h_energy_mc = static_cast<TH1D*>(_file_mc->Get("h_energy"));
    h_energy_mc->SetName("h_energy_mc");
    h_energy_mc->SetDirectory(0);
    h_energy_mc->SetLineWidth(3);
    h_energy_mc->SetLineColor(kRed);
    auto h_energy_corr_mc = static_cast<TH1D*>(_file_mc->Get("h_energy_corr"));
    h_energy_corr_mc->SetName("h_energy_corr_mc");
    h_energy_corr_mc->SetDirectory(0);
    h_energy_corr_mc->SetLineWidth(3);
    h_energy_corr_mc->SetLineColor(kRed);
    auto h_sumRms_mc = static_cast<TH1D*>(_file_mc->Get("h_sumRms"));
    h_sumRms_mc->SetName("h_sumRms_mc");
    h_sumRms_mc->SetDirectory(0);
    h_sumRms_mc->SetLineWidth(3);
    h_sumRms_mc->SetLineColor(kRed);
    auto h_fracLast_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLast"));
    h_fracLast_mc->SetName("h_fracLast_mc");
    h_fracLast_mc->SetDirectory(0);
    h_fracLast_mc->SetLineWidth(3);
    h_fracLast_mc->SetLineColor(kRed);

    auto h_rmsLayer_1_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_1"));
    h_rmsLayer_1_mc->SetName("h_rmsLayer_1_mc");
    h_rmsLayer_1_mc->SetDirectory(0);
    h_rmsLayer_1_mc->SetLineWidth(3);
    h_rmsLayer_1_mc->SetLineColor(kRed);
    auto h_rmsLayer_2_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_2"));
    h_rmsLayer_2_mc->SetName("h_rmsLayer_2_mc");
    h_rmsLayer_2_mc->SetDirectory(0);
    h_rmsLayer_2_mc->SetLineWidth(3);
    h_rmsLayer_2_mc->SetLineColor(kRed);
    auto h_rmsLayer_3_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_3"));
    h_rmsLayer_3_mc->SetName("h_rmsLayer_3_mc");
    h_rmsLayer_3_mc->SetDirectory(0);
    h_rmsLayer_3_mc->SetLineWidth(3);
    h_rmsLayer_3_mc->SetLineColor(kRed);
    auto h_rmsLayer_4_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_4"));
    h_rmsLayer_4_mc->SetName("h_rmsLayer_4_mc");
    h_rmsLayer_4_mc->SetDirectory(0);
    h_rmsLayer_4_mc->SetLineWidth(3);
    h_rmsLayer_4_mc->SetLineColor(kRed);
    auto h_rmsLayer_5_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_5"));
    h_rmsLayer_5_mc->SetName("h_rmsLayer_5_mc");
    h_rmsLayer_5_mc->SetDirectory(0);
    h_rmsLayer_5_mc->SetLineWidth(3);
    h_rmsLayer_5_mc->SetLineColor(kRed);
    auto h_rmsLayer_6_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_6"));
    h_rmsLayer_6_mc->SetName("h_rmsLayer_6_mc");
    h_rmsLayer_6_mc->SetDirectory(0);
    h_rmsLayer_6_mc->SetLineWidth(3);
    h_rmsLayer_6_mc->SetLineColor(kRed);
    auto h_rmsLayer_7_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_7"));
    h_rmsLayer_7_mc->SetName("h_rmsLayer_7_mc");
    h_rmsLayer_7_mc->SetDirectory(0);
    h_rmsLayer_7_mc->SetLineWidth(3);
    h_rmsLayer_7_mc->SetLineColor(kRed);
    auto h_rmsLayer_8_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_8"));
    h_rmsLayer_8_mc->SetName("h_rmsLayer_8_mc");
    h_rmsLayer_8_mc->SetDirectory(0);
    h_rmsLayer_8_mc->SetLineWidth(3);
    h_rmsLayer_8_mc->SetLineColor(kRed);
    auto h_rmsLayer_9_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_9"));
    h_rmsLayer_9_mc->SetName("h_rmsLayer_9_mc");
    h_rmsLayer_9_mc->SetDirectory(0);
    h_rmsLayer_9_mc->SetLineWidth(3);
    h_rmsLayer_9_mc->SetLineColor(kRed);
    auto h_rmsLayer_10_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_10"));
    h_rmsLayer_10_mc->SetName("h_rmsLayer_10_mc");
    h_rmsLayer_10_mc->SetDirectory(0);
    h_rmsLayer_10_mc->SetLineWidth(3);
    h_rmsLayer_10_mc->SetLineColor(kRed);
    auto h_rmsLayer_11_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_11"));
    h_rmsLayer_11_mc->SetName("h_rmsLayer_11_mc");
    h_rmsLayer_11_mc->SetDirectory(0);
    h_rmsLayer_11_mc->SetLineWidth(3);
    h_rmsLayer_11_mc->SetLineColor(kRed);
    auto h_rmsLayer_12_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_12"));
    h_rmsLayer_12_mc->SetName("h_rmsLayer_12_mc");
    h_rmsLayer_12_mc->SetDirectory(0);
    h_rmsLayer_12_mc->SetLineWidth(3);
    h_rmsLayer_12_mc->SetLineColor(kRed);
    auto h_rmsLayer_13_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_13"));
    h_rmsLayer_13_mc->SetName("h_rmsLayer_13_mc");
    h_rmsLayer_13_mc->SetDirectory(0);
    h_rmsLayer_13_mc->SetLineWidth(3);
    h_rmsLayer_13_mc->SetLineColor(kRed);
    auto h_rmsLayer_14_mc = static_cast<TH1D*>(_file_mc->Get("h_rmsLayer_14"));
    h_rmsLayer_14_mc->SetName("h_rmsLayer_14_mc");
    h_rmsLayer_14_mc->SetDirectory(0);
    h_rmsLayer_14_mc->SetLineWidth(3);
    h_rmsLayer_14_mc->SetLineColor(kRed);

    auto h_fracLayer_1_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_1"));
    h_fracLayer_1_mc->SetName("h_fracLayer_1_mc");
    h_fracLayer_1_mc->SetDirectory(0);
    h_fracLayer_1_mc->SetLineWidth(3);
    h_fracLayer_1_mc->SetLineColor(kRed);
    auto h_fracLayer_2_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_2"));
    h_fracLayer_2_mc->SetName("h_fracLayer_2_mc");
    h_fracLayer_2_mc->SetDirectory(0);
    h_fracLayer_2_mc->SetLineWidth(3);
    h_fracLayer_2_mc->SetLineColor(kRed);
    auto h_fracLayer_3_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_3"));
    h_fracLayer_3_mc->SetName("h_fracLayer_3_mc");
    h_fracLayer_3_mc->SetDirectory(0);
    h_fracLayer_3_mc->SetLineWidth(3);
    h_fracLayer_3_mc->SetLineColor(kRed);
    auto h_fracLayer_4_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_4"));
    h_fracLayer_4_mc->SetName("h_fracLayer_4_mc");
    h_fracLayer_4_mc->SetDirectory(0);
    h_fracLayer_4_mc->SetLineWidth(3);
    h_fracLayer_4_mc->SetLineColor(kRed);
    auto h_fracLayer_5_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_5"));
    h_fracLayer_5_mc->SetName("h_fracLayer_5_mc");
    h_fracLayer_5_mc->SetDirectory(0);
    h_fracLayer_5_mc->SetLineWidth(3);
    h_fracLayer_5_mc->SetLineColor(kRed);
    auto h_fracLayer_6_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_6"));
    h_fracLayer_6_mc->SetName("h_fracLayer_6_mc");
    h_fracLayer_6_mc->SetDirectory(0);
    h_fracLayer_6_mc->SetLineWidth(3);
    h_fracLayer_6_mc->SetLineColor(kRed);
    auto h_fracLayer_7_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_7"));
    h_fracLayer_7_mc->SetName("h_fracLayer_7_mc");
    h_fracLayer_7_mc->SetDirectory(0);
    h_fracLayer_7_mc->SetLineWidth(3);
    h_fracLayer_7_mc->SetLineColor(kRed);
    auto h_fracLayer_8_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_8"));
    h_fracLayer_8_mc->SetName("h_fracLayer_8_mc");
    h_fracLayer_8_mc->SetDirectory(0);
    h_fracLayer_8_mc->SetLineWidth(3);
    h_fracLayer_8_mc->SetLineColor(kRed);
    auto h_fracLayer_9_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_9"));
    h_fracLayer_9_mc->SetName("h_fracLayer_9_mc");
    h_fracLayer_9_mc->SetDirectory(0);
    h_fracLayer_9_mc->SetLineWidth(3);
    h_fracLayer_9_mc->SetLineColor(kRed);
    auto h_fracLayer_10_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_10"));
    h_fracLayer_10_mc->SetName("h_fracLayer_10_mc");
    h_fracLayer_10_mc->SetDirectory(0);
    h_fracLayer_10_mc->SetLineWidth(3);
    h_fracLayer_10_mc->SetLineColor(kRed);
    auto h_fracLayer_11_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_11"));
    h_fracLayer_11_mc->SetName("h_fracLayer_11_mc");
    h_fracLayer_11_mc->SetDirectory(0);
    h_fracLayer_11_mc->SetLineWidth(3);
    h_fracLayer_11_mc->SetLineColor(kRed);
    auto h_fracLayer_12_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_12"));
    h_fracLayer_12_mc->SetName("h_fracLayer_12_mc");
    h_fracLayer_12_mc->SetDirectory(0);
    h_fracLayer_12_mc->SetLineWidth(3);
    h_fracLayer_12_mc->SetLineColor(kRed);
    auto h_fracLayer_13_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_13"));
    h_fracLayer_13_mc->SetName("h_fracLayer_13_mc");
    h_fracLayer_13_mc->SetDirectory(0);
    h_fracLayer_13_mc->SetLineWidth(3);
    h_fracLayer_13_mc->SetLineColor(kRed);
    auto h_fracLayer_14_mc = static_cast<TH1D*>(_file_mc->Get("h_fracLayer_14"));
    h_fracLayer_14_mc->SetName("h_fracLayer_14_mc");
    h_fracLayer_14_mc->SetDirectory(0);
    h_fracLayer_14_mc->SetLineWidth(3);
    h_fracLayer_14_mc->SetLineColor(kRed);

    auto h_lastbgolayer_mc = static_cast<TH1D*>(_file_mc->Get("h_lastbgolayer"));
    h_lastbgolayer_mc->SetName("h_lastbgolayer_mc");
    h_lastbgolayer_mc->SetDirectory(0);
    h_lastbgolayer_mc->SetLineWidth(3);
    h_lastbgolayer_mc->SetLineColor(kRed);

    auto h_nbgoentries_mc = static_cast<TH1D*>(_file_mc->Get("h_nbgoentries"));
    h_nbgoentries_mc->SetName("h_nbgoentries_mc");
    h_nbgoentries_mc->SetDirectory(0);
    h_nbgoentries_mc->SetLineWidth(3);
    h_nbgoentries_mc->SetLineColor(kRed);

    auto h_xtrl_mc = static_cast<TH1D*>(_file_mc->Get("h_xtrl"));
    h_xtrl_mc->SetName("h_xtrl_mc");
    h_xtrl_mc->SetDirectory(0);
    h_xtrl_mc->SetLineWidth(3);
    h_xtrl_mc->SetLineColor(kRed);

    auto h_nudadc_1_mc = static_cast<TH1D*>(_file_mc->Get("h_nudadc_1"));
    h_nudadc_1_mc->SetName("h_nudadc_1_mc");
    h_nudadc_1_mc->SetDirectory(0);
    h_nudadc_1_mc->SetLineWidth(3);
    h_nudadc_1_mc->SetLineColor(kRed);
    auto h_nudadc_2_mc = static_cast<TH1D*>(_file_mc->Get("h_nudadc_2"));
    h_nudadc_2_mc->SetName("h_nudadc_2_mc");
    h_nudadc_2_mc->SetDirectory(0);
    h_nudadc_2_mc->SetLineWidth(3);
    h_nudadc_2_mc->SetLineColor(kRed);
    auto h_nudadc_3_mc = static_cast<TH1D*>(_file_mc->Get("h_nudadc_3"));
    h_nudadc_3_mc->SetName("h_nudadc_3_mc");
    h_nudadc_3_mc->SetDirectory(0);
    h_nudadc_3_mc->SetLineWidth(3);
    h_nudadc_3_mc->SetLineColor(kRed);
    auto h_nudadc_4_mc = static_cast<TH1D*>(_file_mc->Get("h_nudadc_4"));
    h_nudadc_4_mc->SetName("h_nudadc_4_mc");
    h_nudadc_4_mc->SetDirectory(0);
    h_nudadc_4_mc->SetLineWidth(3);
    h_nudadc_4_mc->SetLineColor(kRed);

    auto h_nud_totaladc_mc = static_cast<TH1D*>(_file_mc->Get("h_nud_totaladc"));
    h_nud_totaladc_mc->SetName("h_nud_totaladc_mc");
    h_nud_totaladc_mc->SetDirectory(0);
    h_nud_totaladc_mc->SetLineWidth(3);
    h_nud_totaladc_mc->SetLineColor(kRed);
    auto h_nud_maxadc_mc = static_cast<TH1D*>(_file_mc->Get("h_nud_maxadc"));
    h_nud_maxadc_mc->SetName("h_nud_maxadc_mc");
    h_nud_maxadc_mc->SetDirectory(0);
    h_nud_maxadc_mc->SetLineWidth(3);
    h_nud_maxadc_mc->SetLineColor(kRed);

    _file_mc->Close();

    TFile *outfile = TFile::Open(outputfile, "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outputfile << "]" << std::endl;
        exit(100);
    }

    outfile->mkdir("data");
    outfile->cd("data");

    h_energy_data->Write();
    h_energy_corr_data->Write();
    h_sumRms_data->Write();
    h_fracLast_data->Write();
    h_rmsLayer_1_data->Write();
    h_rmsLayer_2_data->Write();
    h_rmsLayer_3_data->Write();
    h_rmsLayer_4_data->Write();
    h_rmsLayer_5_data->Write();
    h_rmsLayer_6_data->Write();
    h_rmsLayer_7_data->Write();
    h_rmsLayer_8_data->Write();
    h_rmsLayer_9_data->Write();
    h_rmsLayer_10_data->Write();
    h_rmsLayer_11_data->Write();
    h_rmsLayer_12_data->Write();
    h_rmsLayer_13_data->Write();
    h_rmsLayer_14_data->Write();
    h_fracLayer_1_data->Write();
    h_fracLayer_2_data->Write();
    h_fracLayer_3_data->Write();
    h_fracLayer_4_data->Write();
    h_fracLayer_5_data->Write();
    h_fracLayer_6_data->Write();
    h_fracLayer_7_data->Write();
    h_fracLayer_8_data->Write();
    h_fracLayer_9_data->Write();
    h_fracLayer_10_data->Write();
    h_fracLayer_11_data->Write();
    h_fracLayer_12_data->Write();
    h_fracLayer_13_data->Write();
    h_fracLayer_14_data->Write();
    h_lastbgolayer_data->Write();
    h_nbgoentries_data->Write();
    h_xtrl_data->Write();
    h_nudadc_1_data->Write();
    h_nudadc_2_data->Write();
    h_nudadc_3_data->Write();
    h_nudadc_4_data->Write();
    h_nud_totaladc_data->Write();
    h_nud_maxadc_data->Write();

    outfile->mkdir("mc");
    outfile->cd("mc");

    h_energy_mc->Write();
    h_energy_corr_mc->Write();
    h_sumRms_mc->Write();
    h_fracLast_mc->Write();
    h_rmsLayer_1_mc->Write();
    h_rmsLayer_2_mc->Write();
    h_rmsLayer_3_mc->Write();
    h_rmsLayer_4_mc->Write();
    h_rmsLayer_5_mc->Write();
    h_rmsLayer_6_mc->Write();
    h_rmsLayer_7_mc->Write();
    h_rmsLayer_8_mc->Write();
    h_rmsLayer_9_mc->Write();
    h_rmsLayer_10_mc->Write();
    h_rmsLayer_11_mc->Write();
    h_rmsLayer_12_mc->Write();
    h_rmsLayer_13_mc->Write();
    h_rmsLayer_14_mc->Write();
    h_fracLayer_1_mc->Write();
    h_fracLayer_2_mc->Write();
    h_fracLayer_3_mc->Write();
    h_fracLayer_4_mc->Write();
    h_fracLayer_5_mc->Write();
    h_fracLayer_6_mc->Write();
    h_fracLayer_7_mc->Write();
    h_fracLayer_8_mc->Write();
    h_fracLayer_9_mc->Write();
    h_fracLayer_10_mc->Write();
    h_fracLayer_11_mc->Write();
    h_fracLayer_12_mc->Write();
    h_fracLayer_13_mc->Write();
    h_fracLayer_14_mc->Write();
    h_lastbgolayer_mc->Write();
    h_nbgoentries_mc->Write();
    h_xtrl_mc->Write();
    h_nudadc_1_mc->Write();
    h_nudadc_2_mc->Write();
    h_nudadc_3_mc->Write();
    h_nudadc_4_mc->Write();
    h_nud_totaladc_mc->Write();
    h_nud_maxadc_mc->Write();

    outfile->mkdir("canvas");
    outfile->cd("canvas");

    TCanvas c_energy("c_energy");
    h_energy_data->DrawNormalized("hist");
    h_energy_mc->DrawNormalized("hist, same");
    c_energy.Write();
    
    TCanvas c_energy_corr("c_energy_corr");
    h_energy_corr_data->DrawNormalized("hist");
    h_energy_corr_mc->DrawNormalized("hist, same");
    c_energy_corr.Write();

    TCanvas c_sumrms("c_sumrms");
    h_sumRms_data->DrawNormalized("hist");
    h_sumRms_mc->DrawNormalized("hist, same");
    c_sumrms.Write();

    TCanvas c_fraclast("c_fraclast");
    h_fracLast_data->DrawNormalized("hist");
    h_fracLast_mc->DrawNormalized("hist, same");
    c_fraclast.Write();

    TCanvas c_rmslayer("c_rmslayer");
    c_rmslayer.Divide(7,2);
    c_rmslayer.cd(1);
    h_rmsLayer_1_data->DrawNormalized("hist");
    h_rmsLayer_1_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(2);
    h_rmsLayer_2_data->DrawNormalized("hist");
    h_rmsLayer_2_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(3);
    h_rmsLayer_3_data->DrawNormalized("hist");
    h_rmsLayer_3_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(4);
    h_rmsLayer_4_data->DrawNormalized("hist");
    h_rmsLayer_4_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(5);
    h_rmsLayer_5_data->DrawNormalized("hist");
    h_rmsLayer_5_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(6);
    h_rmsLayer_6_data->DrawNormalized("hist");
    h_rmsLayer_6_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(7);
    h_rmsLayer_7_data->DrawNormalized("hist");
    h_rmsLayer_7_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(8);
    h_rmsLayer_8_data->DrawNormalized("hist");
    h_rmsLayer_8_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(9);
    h_rmsLayer_9_data->DrawNormalized("hist");
    h_rmsLayer_9_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(10);
    h_rmsLayer_10_data->DrawNormalized("hist");
    h_rmsLayer_10_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(11);
    h_rmsLayer_11_data->DrawNormalized("hist");
    h_rmsLayer_11_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(12);
    h_rmsLayer_12_data->DrawNormalized("hist");
    h_rmsLayer_12_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(13);
    h_rmsLayer_13_data->DrawNormalized("hist");
    h_rmsLayer_13_mc->DrawNormalized("hist, same");

    c_rmslayer.cd(14);
    h_rmsLayer_14_data->DrawNormalized("hist");
    h_rmsLayer_14_mc->DrawNormalized("hist, same");
    c_rmslayer.Write();

    TCanvas c_fraclayer("c_fraclayer");
    c_fraclayer.Divide(7,2);
    c_fraclayer.cd(1);
    h_fracLayer_1_data->DrawNormalized("hist");
    h_fracLayer_1_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(2);
    h_fracLayer_2_data->DrawNormalized("hist");
    h_fracLayer_2_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(3);
    h_fracLayer_3_data->DrawNormalized("hist");
    h_fracLayer_3_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(4);
    h_fracLayer_4_data->DrawNormalized("hist");
    h_fracLayer_4_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(5);
    h_fracLayer_5_data->DrawNormalized("hist");
    h_fracLayer_5_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(6);
    h_fracLayer_6_data->DrawNormalized("hist");
    h_fracLayer_6_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(7);
    h_fracLayer_7_data->DrawNormalized("hist");
    h_fracLayer_7_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(8);
    h_fracLayer_8_data->DrawNormalized("hist");
    h_fracLayer_8_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(9);
    h_fracLayer_9_data->DrawNormalized("hist");
    h_fracLayer_9_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(10);
    h_fracLayer_10_data->DrawNormalized("hist");
    h_fracLayer_10_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(11);
    h_fracLayer_11_data->DrawNormalized("hist");
    h_fracLayer_11_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(12);
    h_fracLayer_12_data->DrawNormalized("hist");
    h_fracLayer_12_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(13);
    h_fracLayer_13_data->DrawNormalized("hist");
    h_fracLayer_13_mc->DrawNormalized("hist, same");

    c_fraclayer.cd(14);
    h_fracLayer_14_data->DrawNormalized("hist");
    h_fracLayer_14_mc->DrawNormalized("hist, same");
    c_fraclayer.Write();

    TCanvas c_lastbgolayer("c_lastbgolayer");
    h_lastbgolayer_data->DrawNormalized("hist");
    h_lastbgolayer_mc->DrawNormalized("hist, same");
    c_lastbgolayer.Write();

    TCanvas c_nbgoentries("c_nbgoentries");
    h_nbgoentries_data->DrawNormalized("hist");
    h_nbgoentries_mc->DrawNormalized("hist, same");
    c_nbgoentries.Write();

    TCanvas c_xtrl("c_xtrl");
    h_xtrl_data->DrawNormalized("hist");
    h_xtrl_mc->DrawNormalized("hist, same");
    c_xtrl.Write();

    TCanvas c_nudadc("c_nudadc");
    c_nudadc.Divide(2,2);
    
    c_nudadc.cd(1);
    h_nudadc_1_data->DrawNormalized("hist");
    h_nudadc_1_mc->DrawNormalized("hist, same");

    c_nudadc.cd(2);
    h_nudadc_2_data->DrawNormalized("hist");
    h_nudadc_2_mc->DrawNormalized("hist, same");

    c_nudadc.cd(3);
    h_nudadc_3_data->DrawNormalized("hist");
    h_nudadc_3_mc->DrawNormalized("hist, same");

    c_nudadc.cd(4);
    h_nudadc_4_data->DrawNormalized("hist");
    h_nudadc_4_mc->DrawNormalized("hist, same");
    c_nudadc.Write();

    TCanvas c_totaladc("c_totaladc");
    h_nud_totaladc_data->DrawNormalized("hist");
    h_nud_totaladc_mc->DrawNormalized("hist, same");
    c_totaladc.Write();

    TCanvas c_maxadc("c_maxadc");
    h_nud_maxadc_data->DrawNormalized("hist");
    h_nud_maxadc_mc->DrawNormalized("hist, same");
    c_maxadc.Write();

    outfile->Close();
}