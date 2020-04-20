#include "myHeader.h"
#include "acceptance.h"

void generateFinalGraph(
    const bool verbose,
    const bool pedantic,
    const std::string outputPath,
    const std::string complete_histo_path,
    const std::string wd)
{
    // Read energy log-binning
    auto logEBins = readLogBinning(wd);
    if (pedantic)
    {
        std::cout << "\nEnergy log binning..." << std::scientific;
        for (auto it = logEBins.begin(); it != logEBins.end(); ++it)
            std::cout << "\n"
                      << *it;
        std::cout << std::defaultfloat;
    }

    // Create and load acceptance events cuts from config file
    acceptance_conf acceptance_cuts;
    acceptance_active_cuts active_cuts;
    load_acceptance_struct(acceptance_cuts, active_cuts, wd);

    // Open input complete histos
    TFile inHisto(complete_histo_path.c_str(), "READ");
    if (!inHisto.IsOpen())
    {
        std::cerr << "Error reading input complete histo: " << complete_histo_path;
        exit(123);
    }
    
    auto h_incoming = static_cast<TH1D *>(inHisto.Get("h_incoming"));
    auto h_gometric_cut = static_cast<TH1D *>(inHisto.Get("h_gometric_cut"));
    auto h_maxElayer_cut = static_cast<TH1D *>(inHisto.Get("h_maxElayer_cut"));
    auto h_maxBarLayer_cut = static_cast<TH1D *>(inHisto.Get("h_maxBarLayer_cut"));
    auto h_BGOTrackContainment_cut = static_cast<TH1D *>(inHisto.Get("h_BGOTrackContainment_cut"));
    auto h_BGO_fiducial = static_cast<TH1D *>(inHisto.Get("h_BGO_fiducial_cut"));
    auto h_nBarLayer13_cut = static_cast<TH1D *>(inHisto.Get("h_nBarLayer13_cut"));
    auto h_maxRms_cut = static_cast<TH1D *>(inHisto.Get("h_maxRms_cut"));
    auto h_track_selection_cut = static_cast<TH1D *>(inHisto.Get("h_track_selection_cut"));
    auto h_xtrl_cut = static_cast<TH1D *>(inHisto.Get("h_xtrl_cut"));
    auto h_psd_charge_cut = static_cast<TH1D *>(inHisto.Get("h_psd_charge_cut"));
    auto h_all_cut = static_cast<TH1D *>(inHisto.Get("h_all_cut"));

    h_incoming->SetDirectory(0);
    h_gometric_cut->SetDirectory(0);
    h_maxElayer_cut->SetDirectory(0);
    h_maxBarLayer_cut->SetDirectory(0);
    h_BGOTrackContainment_cut->SetDirectory(0);
    h_BGO_fiducial->SetDirectory(0);
    h_nBarLayer13_cut->SetDirectory(0);
    h_maxRms_cut->SetDirectory(0);
    h_track_selection_cut->SetDirectory(0);
    h_xtrl_cut->SetDirectory(0);
    h_psd_charge_cut->SetDirectory(0);
    h_all_cut->SetDirectory(0);

    inHisto.Close();
    
    // Building acceptance histos
    auto h_acceptance_gometric_cut = static_cast<TH1D *>(h_gometric_cut->Clone("h_acceptance_gometric_cut"));
    auto h_acceptance_maxElayer_cut = static_cast<TH1D *>(h_maxElayer_cut->Clone("h_acceptance_maxElayer_cut"));
    auto h_acceptance_maxBarLayer_cut = static_cast<TH1D *>(h_maxBarLayer_cut->Clone("h_acceptance_maxBarLayer_cut"));
    auto h_acceptance_BGOTrackContainment_cut = static_cast<TH1D *>(h_BGOTrackContainment_cut->Clone("h_acceptance_BGOTrackContainment_cut"));
    auto h_acceptance_BGO_fiducial = static_cast<TH1D *>(h_BGO_fiducial->Clone("h_acceptance_BGO_fiducial"));
    auto h_acceptance_nBarLayer13_cut = static_cast<TH1D *>(h_nBarLayer13_cut->Clone("h_acceptance_nBarLayer13_cut"));
    auto h_acceptance_maxRms_cut = static_cast<TH1D *>(h_maxRms_cut->Clone("h_acceptance_maxRms_cut"));
    auto h_acceptance_track_selection_cut = static_cast<TH1D *>(h_track_selection_cut->Clone("h_acceptance_track_selection_cut"));
    auto h_acceptance_xtrl_cut = static_cast<TH1D *>(h_xtrl_cut->Clone("h_acceptance_xtrl_cut"));
    auto h_acceptance_psd_charge_cut = static_cast<TH1D *>(h_psd_charge_cut->Clone("h_acceptance_psd_charge_cut"));
    auto h_acceptance_all_cut = static_cast<TH1D *>(h_all_cut->Clone("h_acceptance_all_cut"));

    // Building ratio histos
    auto h_ratio_gometric_cut = static_cast<TH1D *>(h_gometric_cut->Clone("h_ratio_gometric_cut"));
    auto h_ratio_maxElayer_cut = static_cast<TH1D *>(h_maxElayer_cut->Clone("h_ratio_maxElayer_cut"));
    auto h_ratio_maxBarLayer_cut = static_cast<TH1D *>(h_maxBarLayer_cut->Clone("h_ratio_maxBarLayer_cut"));
    auto h_ratio_BGOTrackContainment_cut = static_cast<TH1D *>(h_BGOTrackContainment_cut->Clone("h_ratio_BGOTrackContainment_cut"));
    auto h_ratio_BGO_fiducial = static_cast<TH1D *>(h_BGO_fiducial->Clone("h_ratio_BGO_fiducial"));
    auto h_ratio_nBarLayer13_cut = static_cast<TH1D *>(h_nBarLayer13_cut->Clone("h_ratio_nBarLayer13_cut"));
    auto h_ratio_maxRms_cut = static_cast<TH1D *>(h_maxRms_cut->Clone("h_ratio_maxRms_cut"));
    auto h_ratio_track_selection_cut = static_cast<TH1D *>(h_track_selection_cut->Clone("h_ratio_track_selection_cut"));
    auto h_ratio_xtrl_cut = static_cast<TH1D *>(h_xtrl_cut->Clone("h_ratio_xtrl_cut"));
    auto h_ratio_psd_charge_cut = static_cast<TH1D *>(h_psd_charge_cut->Clone("h_ratio_psd_charge_cut"));
    auto h_ratio_all_cut = static_cast<TH1D *>(h_all_cut->Clone("h_ratio_all_cut"));
   
    double genSurface = 4 * TMath::Pi() * pow(acceptance_cuts.vertex_radius, 2) / 2;

    h_acceptance_gometric_cut->Scale(genSurface);
    h_acceptance_maxElayer_cut->Scale(genSurface);
    h_acceptance_maxBarLayer_cut->Scale(genSurface);
    h_acceptance_BGOTrackContainment_cut->Scale(genSurface);
    h_acceptance_BGO_fiducial->Scale(genSurface);
    h_acceptance_nBarLayer13_cut->Scale(genSurface);
    h_acceptance_maxRms_cut->Scale(genSurface);
    h_acceptance_track_selection_cut->Scale(genSurface);
    h_acceptance_xtrl_cut->Scale(genSurface);
    h_acceptance_psd_charge_cut->Scale(genSurface);
    h_acceptance_all_cut->Scale(genSurface);

    h_acceptance_gometric_cut->Divide(h_incoming);
    h_acceptance_maxElayer_cut->Divide(h_incoming);
    h_acceptance_maxBarLayer_cut->Divide(h_incoming);
    h_acceptance_BGOTrackContainment_cut->Divide(h_incoming);
    h_acceptance_BGO_fiducial->Divide(h_incoming);
    h_acceptance_nBarLayer13_cut->Divide(h_incoming);
    h_acceptance_maxRms_cut->Divide(h_incoming);
    h_acceptance_track_selection_cut->Divide(h_incoming);
    h_acceptance_xtrl_cut->Divide(h_incoming);
    h_acceptance_psd_charge_cut->Divide(h_incoming);
    h_acceptance_all_cut->Divide(h_incoming);

    h_ratio_gometric_cut->Divide(h_incoming);
    h_ratio_maxElayer_cut->Divide(h_incoming);
    h_ratio_maxBarLayer_cut->Divide(h_incoming);
    h_ratio_BGOTrackContainment_cut->Divide(h_incoming);
    h_ratio_BGO_fiducial->Divide(h_incoming);
    h_ratio_nBarLayer13_cut->Divide(h_incoming);
    h_ratio_maxRms_cut->Divide(h_incoming);
    h_ratio_track_selection_cut->Divide(h_incoming);
    h_ratio_xtrl_cut->Divide(h_incoming);
    h_ratio_psd_charge_cut->Divide(h_incoming);
    h_ratio_all_cut->Divide(h_incoming);

    // Builing vectors
    std::vector<double> energyValues(h_incoming->GetXaxis()->GetNbins(), 0);

    std::vector<double> acceptanceValues_gometric_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxElayer_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxBarLayer_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGOTrackContainment_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGO_fiducial_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_nBarLayer13_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxRms_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_track_selection_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_xtrl_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_psd_charge_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_all_cut(energyValues.size(), 0);

    for (auto it = logEBins.begin(); it != (logEBins.end() - 1); ++it)
    {
        auto index = std::distance(logEBins.begin(), it);
        energyValues[index] = wtsydp(*it, *(it + 1), -1);
        acceptanceValues_gometric_cut[index] = h_acceptance_gometric_cut->GetBinContent(index + 1);
        acceptanceValues_maxElayer_cut[index] = h_acceptance_maxElayer_cut->GetBinContent(index + 1);
        acceptanceValues_maxBarLayer_cut[index] = h_acceptance_maxBarLayer_cut->GetBinContent(index + 1);
        acceptanceValues_BGOTrackContainment_cut[index] = h_acceptance_BGOTrackContainment_cut->GetBinContent(index + 1);
        acceptanceValues_BGO_fiducial_cut[index] = h_acceptance_BGO_fiducial->GetBinContent(index + 1);
        acceptanceValues_nBarLayer13_cut[index] = h_acceptance_nBarLayer13_cut->GetBinContent(index + 1);
        acceptanceValues_maxRms_cut[index] = h_acceptance_maxRms_cut->GetBinContent(index + 1);
        acceptanceValues_track_selection_cut[index] = h_acceptance_track_selection_cut->GetBinContent(index + 1);
        acceptanceValues_xtrl_cut[index] = h_acceptance_xtrl_cut->GetBinContent(index + 1);
        acceptanceValues_psd_charge_cut[index] = h_acceptance_psd_charge_cut->GetBinContent(index + 1);
        acceptanceValues_all_cut[index] = h_acceptance_all_cut->GetBinContent(index + 1);
    }

    // Building graphs
    TGraph gr_acceptance_gometric_cut(energyValues.size(), &energyValues[0], &acceptanceValues_gometric_cut[0]);
    TGraph gr_acceptance_maxElayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxElayer_cut[0]);
    TGraph gr_acceptance_maxBarLayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxBarLayer_cut[0]);
    TGraph gr_acceptance_BGOTrackContainment_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGOTrackContainment_cut[0]);
    TGraph gr_acceptance_BGO_fiducial_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_cut[0]);
    TGraph gr_acceptance_nBarLayer13_cut(energyValues.size(), &energyValues[0], &acceptanceValues_nBarLayer13_cut[0]);
    TGraph gr_acceptance_maxRms_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxRms_cut[0]);
    TGraph gr_acceptance_track_selection_cut(energyValues.size(), &energyValues[0], &acceptanceValues_track_selection_cut[0]);
    TGraph gr_acceptance_xtrl_cut(energyValues.size(), &energyValues[0], &acceptanceValues_xtrl_cut[0]);
    TGraph gr_acceptance_psd_charge_cut(energyValues.size(), &energyValues[0], &acceptanceValues_psd_charge_cut[0]);
    TGraph gr_acceptance_all_cut(energyValues.size(), &energyValues[0], &acceptanceValues_all_cut[0]);

    gr_acceptance_gometric_cut.SetName("gr_acceptance_gometric_cut");
    gr_acceptance_maxElayer_cut.SetName("gr_acceptance_maxElayer_cut");
    gr_acceptance_maxBarLayer_cut.SetName("gr_acceptance_maxBarLayer_cut");
    gr_acceptance_BGOTrackContainment_cut.SetName("gr_acceptance_BGOTrackContainment_cut");
    gr_acceptance_BGO_fiducial_cut.SetName("gr_acceptance_BGO_fiducial_cut");
    gr_acceptance_nBarLayer13_cut.SetName("gr_acceptance_nBarLayer13_cut");
    gr_acceptance_maxRms_cut.SetName("gr_acceptance_maxRms_cut");
    gr_acceptance_track_selection_cut.SetName("gr_acceptance_track_selection_cut");
    gr_acceptance_xtrl_cut.SetName("gr_acceptance_xtrl_cut");
    gr_acceptance_psd_charge_cut.SetName("gr_acceptance_psd_charge_cut");
    gr_acceptance_all_cut.SetName("gr_acceptance_all_cut");

    gr_acceptance_gometric_cut.SetTitle("Acceptance - geometric cut");
    gr_acceptance_maxElayer_cut.SetTitle("Acceptance - maxElateral cut");
    gr_acceptance_maxBarLayer_cut.SetTitle("Acceptance - maxBarLayer cut");
    gr_acceptance_BGOTrackContainment_cut.SetTitle("Acceptance - BGOTrackContainment cut");
    gr_acceptance_BGO_fiducial_cut.SetTitle("Acceptance - BGO fiducial volume cut");
    gr_acceptance_nBarLayer13_cut.SetTitle("Acceptance - nBarLayer13 cut");
    gr_acceptance_maxRms_cut.SetTitle("Acceptance - maxRms cut");
    gr_acceptance_track_selection_cut.SetTitle("Acceptance - track selection cut");
    gr_acceptance_xtrl_cut.SetTitle("Acceptance - XTRL cut");
    gr_acceptance_psd_charge_cut.SetTitle("Acceptance - PSD charge selection cut");
    gr_acceptance_all_cut.SetTitle("Acceptance - all cut");

    // Write output TFile
    TFile outFile(outputPath.c_str(), "RECREATE");
    if (!outFile.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outputPath << std::endl;
        exit(123);
    }

    // Write final TGraphs
    gr_acceptance_gometric_cut.Write();
    gr_acceptance_maxElayer_cut.Write();
    gr_acceptance_maxBarLayer_cut.Write();
    gr_acceptance_BGOTrackContainment_cut.Write();
    gr_acceptance_BGO_fiducial_cut.Write();
    gr_acceptance_nBarLayer13_cut.Write();
    gr_acceptance_maxRms_cut.Write();
    gr_acceptance_track_selection_cut.Write();
    gr_acceptance_xtrl_cut.Write();
    gr_acceptance_psd_charge_cut.Write();
    gr_acceptance_all_cut.Write();

    // Write ratio histos
    h_ratio_gometric_cut->Write();
    h_ratio_maxElayer_cut->Write();
    h_ratio_maxBarLayer_cut->Write();
    h_ratio_BGOTrackContainment_cut->Write();
    h_ratio_BGO_fiducial->Write();
    h_ratio_nBarLayer13_cut->Write();
    h_ratio_maxRms_cut->Write();
    h_ratio_track_selection_cut->Write();
    h_ratio_xtrl_cut->Write();
    h_ratio_psd_charge_cut->Write();
    h_ratio_all_cut->Write();

    // Close output TFile
    outFile.Close();
}