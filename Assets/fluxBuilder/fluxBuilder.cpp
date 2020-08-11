#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"

#include <iostream>

std::vector<float> createLogBinning(
	const double eMin,
	const double eMax,
	const std::size_t n_bins);

double wtsydp(
    const float minene,
    const float maxene,
    const float index);

void fluxBuilder(
    const char* electronMC, 
    const char* protonMC, 
    const char* data, 
    const double livetime, 
    const double emin = 1, 
    const double emax = 1e+4)
{   
    // Extract electron acceptance
    TFile* electron_mc_file = TFile::Open(electronMC, "READ");
    if (electron_mc_file->IsZombie())
    {
        std::cerr << "\n\nError reading input electron MC reco file: " << electronMC << std::endl;
        exit(100);
    }

    TH1D* electron_acceptance = static_cast<TH1D*>(electron_mc_file->Get("Acceptance_histos/h_acceptance_all_cut"));
    electron_acceptance->SetDirectory(0);

    electron_mc_file->Close();

    // Extract proton contamination
    TFile* proton_mc_file = TFile::Open(protonMC, "READ");
    if (proton_mc_file->IsZombie())
    {
        std::cerr << "\n\nError reading input proton MC reco file: " << protonMC << std::endl;
        exit(100);
    }

    TH1D* proton_background_fraction = static_cast<TH1D*>(proton_mc_file->Get("mc_ancillary/proton_background_ratio"));
    proton_background_fraction->SetDirectory(0);

    proton_mc_file->Close();

    // Read data file
    TFile* data_file = TFile::Open(data, "READ");
    if (data_file->IsZombie())
    {
        std::cerr << "\n\nError reading input data reco file: " << data << std::endl;
        exit(100);
    }

    TH1D* h_electron_data_all_cut = static_cast<TH1D*>(data_file->Get("h_all_cut"));
    TH1D* h_electron_data_over_xtrl_cut = static_cast<TH1D*>(data_file->Get("mc_ancillary/h_background_over_xtrl_cut"));

    h_electron_data_all_cut->SetDirectory(0);
    h_electron_data_over_xtrl_cut->SetDirectory(0);

    data_file->Close();

    // Compute background protons
    TH1D* backgroung_protons = static_cast<TH1D*>(proton_background_fraction->Clone("backgroung_protons"));
    backgroung_protons->Multiply(h_electron_data_over_xtrl_cut);

    // Compute dN/dE
    TH1D* dNdE = static_cast<TH1D*>(h_electron_data_all_cut->Clone("dNdE"));
    dNdE->Add(backgroung_protons, -1);

    // Rescale by energy bin width
    for (int bIdx=1; bIdx <= dNdE->GetNbinsX(); ++bIdx)
    {
        dNdE->SetBinContent(bIdx, dNdE->GetBinContent(bIdx)/dNdE->GetBinWidth(bIdx));
        dNdE->SetBinError(bIdx, dNdE->GetBinError(bIdx)/dNdE->GetBinWidth(bIdx));
    }

    // Build flux
    TH1D* all_electron_flux = static_cast<TH1D*>(dNdE->Clone("all_electron_flux"));
    all_electron_flux->Divide(electron_acceptance);
    all_electron_flux->Scale(1/livetime);
    
    // Build energy binning
    auto logEBins = createLogBinning(emin, emax, all_electron_flux->GetNbinsX());

    // Build energy vector
    std::vector<double> energyValues (logEBins.size(), 0);
    for (auto it = logEBins.begin(); it != (logEBins.end() - 1); ++it)
        energyValues[std::distance(logEBins.begin(), it)] = wtsydp(*it, *(it + 1), -3);

    // Build flux vectors
    std::vector<double> energy_errors(all_electron_flux->GetNbinsX(), 0);
    std::vector<double> all_electron_flux_values (all_electron_flux->GetNbinsX(), 0);
    std::vector<double> all_electron_flux_errors (all_electron_flux->GetNbinsX(), 0);
    std::vector<double> all_electron_fluxE3_values (all_electron_flux->GetNbinsX(), 0);
    std::vector<double> all_electron_fluxE3_errors (all_electron_flux->GetNbinsX(), 0);

    for (int bIdx=1; bIdx<=all_electron_flux->GetNbinsX(); ++bIdx)
    {
        all_electron_flux_values[bIdx-1] = all_electron_flux->GetBinContent(bIdx);
        all_electron_flux_errors[bIdx-1] = all_electron_flux->GetBinError(bIdx);
    }

    // Build E^3 flux
    TH1D* all_electron_flux_E3 = static_cast<TH1D*>(all_electron_flux->Clone("all_electron_flux_E3"));
    for (int bIdx=1; bIdx<=all_electron_flux_E3->GetNbinsX(); ++bIdx)
    {
        all_electron_flux_E3->SetBinContent(bIdx, all_electron_flux_E3->GetBinContent(bIdx) * pow(all_electron_flux_E3->GetBinCenter(bIdx), 3 )); 
        all_electron_flux_E3->SetBinError(bIdx, all_electron_flux_E3->GetBinError(bIdx) * pow(all_electron_flux_E3->GetBinCenter(bIdx), 3 ));
        all_electron_fluxE3_values[bIdx-1] = all_electron_flux_E3->GetBinContent(bIdx);
        all_electron_fluxE3_errors[bIdx-1] = all_electron_flux_E3->GetBinError(bIdx);
    }

    // Build flux TGraphErrors
    TGraphErrors gr_all_electron_flux(energyValues.size(), &(energyValues[0]), &(all_electron_flux_values[0]), &(energy_errors[0]), &(all_electron_flux_errors[0]));
    TGraphErrors gr_all_electron_flux_E3(energyValues.size(), &(energyValues[0]), &(all_electron_fluxE3_values[0]), &(energy_errors[0]), &(all_electron_fluxE3_errors[0]));

    gr_all_electron_flux.SetName("gr_all_electron_flux");
    gr_all_electron_flux_E3.SetName("gr_all_electron_flux_E3");
    gr_all_electron_flux.SetTitle("All Electron Flux");
    gr_all_electron_flux_E3.SetTitle("All Electron Flux");

    TFile* all_electron_flux_file = TFile::Open("all_electron_flux.root", "RECREATE");
    if (all_electron_flux_file->IsZombie())
    {
        std::cerr << "\n\nError writing final flux ROOT file" << std::endl;
        exit(100);
    }

    electron_acceptance->Write();
    proton_background_fraction->Write();
    h_electron_data_all_cut->Write();
    h_electron_data_over_xtrl_cut->Write();
    backgroung_protons->Write();
    dNdE->Write();
    all_electron_flux->Write();
    all_electron_flux_E3->Write();
    gr_all_electron_flux.Write();
    gr_all_electron_flux_E3.Write();

    all_electron_flux_file->Close();
}

std::vector<float> createLogBinning(
	const double eMin,
	const double eMax,
	const std::size_t n_bins)
{
	std::vector<float> binning(n_bins + 1, 0);
	double log_interval = (log10(eMax) - log10(eMin)) / n_bins;
	for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
		binning[bIdx] = pow(10, log10(eMin) + bIdx * log_interval);

	return binning;
}

double wtsydp(
    const float minene,
    const float maxene,
    const float index)
{
    float dene = maxene - minene;
    if (index != -1)
        return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * dene)), 1. / index);
    else
        return dene / log(maxene / minene);
}