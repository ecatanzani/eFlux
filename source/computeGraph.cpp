#include "myHeader.h"
#include "acceptance.h"

void generateFinalGraph(
    const bool verbose,
    const bool pedantic,
    const std::string outputPath,
    const std::string complete_histo_path,
    const std::string wd)
{
    // Create energy log-binning
    energy_cuts eCuts;
    load_energy_struct(eCuts, wd);
    auto logEBins = createLogBinning(eCuts);
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
    TFile inHisto(complete_histo_path.c_str(),"READ");
    if (!inHisto.IsOpen())
    {
        std::cerr << "Error reading input complete histo: " << complete_histo_path;
        exit(123);
    }
    auto h_incoming = static_cast<TH1D*> (inHisto.Get("h_incoming"));
    auto h_maxElateral_cut = static_cast<TH1D*> (inHisto.Get("h_maxElateral_cut"));
    auto h_maxBarLayer_cut = static_cast<TH1D*> (inHisto.Get("h_maxBarLayer_cut"));
    auto h_BGOTrackContainment_cut = static_cast<TH1D*> (inHisto.Get("h_BGOTrackContainment_cut"));
    auto h_nBarLayer13_cut = static_cast<TH1D*> (inHisto.Get("h_nBarLayer13_cut"));
    auto h_maxRms_cut = static_cast<TH1D*> (inHisto.Get("h_maxRms_cut"));
    auto h_all_cut = static_cast<TH1D*> (inHisto.Get("h_all_cut"));

    h_incoming->SetDirectory(0);
    h_maxElateral_cut->SetDirectory(0);
    h_maxBarLayer_cut->SetDirectory(0);
    h_BGOTrackContainment_cut->SetDirectory(0);
    h_nBarLayer13_cut->SetDirectory(0);
    h_maxRms_cut->SetDirectory(0);
    h_all_cut->SetDirectory(0);
    
    inHisto.Close();

    // Building acceptance histos
    auto h_acceptance_maxElateral_cut = static_cast<TH1D *>(h_maxElateral_cut->Clone("h_acceptance_maxElateral_cut"));
    auto h_acceptance_maxBarLayer_cut = static_cast<TH1D *>(h_maxBarLayer_cut->Clone("h_acceptance_maxBarLayer_cut"));
    auto h_acceptance_BGOTrackContainment_cut = static_cast<TH1D *>(h_BGOTrackContainment_cut->Clone("h_acceptance_BGOTrackContainment_cut"));
    auto h_acceptance_nBarLayer13_cut = static_cast<TH1D *>(h_nBarLayer13_cut->Clone("h_acceptance_nBarLayer13_cut"));
    auto h_acceptance_maxRms_cut = static_cast<TH1D *>(h_maxRms_cut->Clone("h_acceptance_maxRms_cut"));
    auto h_acceptance_all_cut = static_cast<TH1D *>(h_all_cut->Clone("h_acceptance_all_cut"));

    // Generate acceptance histos
    double genSurface = 4 * TMath::Pi() * pow(acceptance_cuts.vertex_radius, 2) / 2;

    h_acceptance_maxElateral_cut->Scale(genSurface);
    h_acceptance_maxBarLayer_cut->Scale(genSurface);
    h_acceptance_BGOTrackContainment_cut->Scale(genSurface);
    h_acceptance_nBarLayer13_cut->Scale(genSurface);
    h_acceptance_maxRms_cut->Scale(genSurface);
    h_acceptance_all_cut->Scale(genSurface);

    h_acceptance_maxElateral_cut->Divide(h_incoming);
    h_acceptance_maxBarLayer_cut->Divide(h_incoming);
    h_acceptance_BGOTrackContainment_cut->Divide(h_incoming);
    h_acceptance_nBarLayer13_cut->Divide(h_incoming);
    h_acceptance_maxRms_cut->Divide(h_incoming);
    h_acceptance_all_cut->Divide(h_incoming);

    // Builing vectors
    std::vector<double> energyValues(h_acceptance_maxElateral_cut->GetXaxis()->GetNbins(), 0);

    std::vector<double> acceptanceValues_maxElateral_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxBarLayer_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_BGOTrackContainment_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_nBarLayer13_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_maxRms_cut(energyValues.size(), 0);
    std::vector<double> acceptanceValues_all_cut(energyValues.size(), 0);

    for (auto it = logEBins.begin(); it != logEBins.end(); ++it)
    {
        auto index = std::distance(logEBins.begin(), it);
        energyValues[index] = wtsydp(*it, *(it + 1), -1);
        acceptanceValues_maxElateral_cut[index] = h_acceptance_maxElateral_cut->GetBinContent(index + 1);
        acceptanceValues_maxBarLayer_cut[index] = h_acceptance_maxBarLayer_cut->GetBinContent(index + 1);
        acceptanceValues_BGOTrackContainment_cut[index] = h_acceptance_BGOTrackContainment_cut->GetBinContent(index + 1);
        acceptanceValues_nBarLayer13_cut[index] = h_acceptance_nBarLayer13_cut->GetBinContent(index + 1);
        acceptanceValues_maxRms_cut[index] = h_acceptance_maxRms_cut->GetBinContent(index + 1);
        acceptanceValues_all_cut[index] = h_acceptance_all_cut->GetBinContent(index + 1);
    }

    // Building graphs
    TGraph gr_acceptance_maxElateral_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxElateral_cut[0]);
    TGraph gr_acceptance_maxBarLayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxBarLayer_cut[0]);
    TGraph gr_acceptance_BGOTrackContainment_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGOTrackContainment_cut[0]);
    TGraph gr_acceptance_nBarLayer13_cut(energyValues.size(), &energyValues[0], &acceptanceValues_nBarLayer13_cut[0]);
    TGraph gr_acceptance_maxRms_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxRms_cut[0]);
    TGraph gr_acceptance_all_cut(energyValues.size(), &energyValues[0], &acceptanceValues_all_cut[0]);

    gr_acceptance_maxElateral_cut.SetName("gr_acceptance_maxElateral_cut");
    gr_acceptance_maxBarLayer_cut.SetName("gr_acceptance_maxBarLayer_cut");
    gr_acceptance_BGOTrackContainment_cut.SetName("gr_acceptance_BGOTrackContainment_cut");
    gr_acceptance_nBarLayer13_cut.SetName("gr_acceptance_nBarLayer13_cut");
    gr_acceptance_maxRms_cut.SetName("gr_acceptance_maxRms_cut");
    gr_acceptance_all_cut.SetName("gr_acceptance_all_cut");

    gr_acceptance_maxElateral_cut.SetTitle("Acceptance - maxElateral cut");
    gr_acceptance_maxBarLayer_cut.SetTitle("Acceptance - maxBarLayer cut");
    gr_acceptance_BGOTrackContainment_cut.SetTitle("Acceptance - BGOTrackContainment cut");
    gr_acceptance_nBarLayer13_cut.SetTitle("Acceptance - nBarLayer13 cut");
    gr_acceptance_maxRms_cut.SetTitle("Acceptance - maxRms cut");
    gr_acceptance_all_cut.SetTitle("Acceptance - all cut");

    // Write output TFile
    TFile outFile(outputPath.c_str(),"RECREATE");
    if (!outFile.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outputPath << std::endl;
        exit(123);
    }

    gr_acceptance_maxElateral_cut.Write();
    gr_acceptance_maxBarLayer_cut.Write();
    gr_acceptance_BGOTrackContainment_cut.Write();
    gr_acceptance_nBarLayer13_cut.Write();
    gr_acceptance_maxRms_cut.Write();
    gr_acceptance_all_cut.Write();

    // Close output TFile
    outFile.Close();
}