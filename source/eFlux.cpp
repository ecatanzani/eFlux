#include "myHeader.h"
#include "acceptance.h"
#include "flux.h"
#include "data_loop.h"
#include "binning.h"
#include "wtsydp.h"

#include "TGraphErrors.h"

void buildFlux(
	const std::string inputPath,
	const unsigned int lvTime,
	TFile &outFile,
	const bool verbose,
	const std::string accInputPath,
	const TH1D *electron_acceptance,
	const TH1D *proton_background_fraction,
	const std::string wd)
{
	/*
    auto acceptance_tf1 = readAcceptance(
        outFile,
        verbose,
        accInputPath);
    */

	// Electron data selection
	auto data_selection =
		evLoop(
			inputPath,
			outFile,
			verbose,
			wd);

#if 0	
	std::unique_ptr<TH1D> h_electron_data_all_cut = std::make_unique<TH1D>(data_selection[0]);
	std::unique_ptr<TH1D> h_electron_data_over_xtrl_cut = std::make_unique<TH1D>(data_selection[1]);

	// Compute background protons
	if (verbose)
		std::cout << "\nComputing proton background...";
	TH1D *backgroung_protons = static_cast<TH1D *>(proton_background_fraction->Clone("backgroung_protons"));
	backgroung_protons->Multiply(h_electron_data_over_xtrl_cut.get());

	// Compute dN/dE
	if (verbose)
		std::cout << "\n\nComputing dN/dE...";
	TH1D *dNdE = static_cast<TH1D *>(h_electron_data_all_cut->Clone("dNdE"));
	dNdE->Add(backgroung_protons, -1);

	// Rescale by energy bin width
	if (verbose)
		std::cout << "\n\nRescaling by energy bin width...";
	for (int bIdx = 1; bIdx <= dNdE->GetNbinsX(); ++bIdx)
	{
		dNdE->SetBinContent(bIdx, dNdE->GetBinContent(bIdx) / dNdE->GetBinWidth(bIdx));
		dNdE->SetBinError(bIdx, dNdE->GetBinError(bIdx) / dNdE->GetBinWidth(bIdx));
	}

	// Build flux
	if (verbose)
		std::cout << "\n\nBuilding flux...";
	TH1D *all_electron_flux = static_cast<TH1D *>(dNdE->Clone("all_electron_flux"));
	all_electron_flux->Divide(electron_acceptance);
	all_electron_flux->Scale(1 / lvTime);

	// Build energy binning
	if (verbose)
		std::cout << "\n\nBuilding energy binning...";
	auto logEBins = createLogBinning(
		all_electron_flux->GetXaxis()->GetXmin(),
		all_electron_flux->GetXaxis()->GetXmax(),
		all_electron_flux->GetNbinsX());

	// Build energy vector
	if (verbose)
		std::cout << "\n\nBuilding energy vectors...";
	std::vector<double> energyValues(logEBins.size(), 0);
	for (auto it = logEBins.begin(); it != (logEBins.end() - 1); ++it)
		energyValues[std::distance(logEBins.begin(), it)] = wtsydp(*it, *(it + 1), -3);

	// Build flux vectors
	if (verbose)
		std::cout << "\n\nBuilding flux vectors...";
	std::vector<double> energy_errors(all_electron_flux->GetNbinsX(), 0);
	std::vector<double> all_electron_flux_values(all_electron_flux->GetNbinsX(), 0);
	std::vector<double> all_electron_flux_errors(all_electron_flux->GetNbinsX(), 0);
	std::vector<double> all_electron_fluxE3_values(all_electron_flux->GetNbinsX(), 0);
	std::vector<double> all_electron_fluxE3_errors(all_electron_flux->GetNbinsX(), 0);

	for (int bIdx = 1; bIdx <= all_electron_flux->GetNbinsX(); ++bIdx)
	{
		all_electron_flux_values[bIdx - 1] = all_electron_flux->GetBinContent(bIdx);
		all_electron_flux_errors[bIdx - 1] = all_electron_flux->GetBinError(bIdx);
	}

	// Build E^3 flux
	if (verbose)
		std::cout << "\n\nBuilding E3 flux...";
	TH1D *all_electron_flux_E3 = static_cast<TH1D *>(all_electron_flux->Clone("all_electron_flux_E3"));
	for (int bIdx = 1; bIdx <= all_electron_flux_E3->GetNbinsX(); ++bIdx)
	{
		all_electron_flux_E3->SetBinContent(bIdx, all_electron_flux_E3->GetBinContent(bIdx) * pow(all_electron_flux_E3->GetBinCenter(bIdx), 3));
		all_electron_flux_E3->SetBinError(bIdx, all_electron_flux_E3->GetBinError(bIdx) * pow(all_electron_flux_E3->GetBinCenter(bIdx), 3));
		all_electron_fluxE3_values[bIdx - 1] = all_electron_flux_E3->GetBinContent(bIdx);
		all_electron_fluxE3_errors[bIdx - 1] = all_electron_flux_E3->GetBinError(bIdx);
	}

	// Build flux TGraphErrors
	if (verbose)
		std::cout << "\n\nBuilding TGraphErrors...";
	TGraphErrors gr_all_electron_flux(energyValues.size(), &(energyValues[0]), &(all_electron_flux_values[0]), &(energy_errors[0]), &(all_electron_flux_errors[0]));
	TGraphErrors gr_all_electron_flux_E3(energyValues.size(), &(energyValues[0]), &(all_electron_fluxE3_values[0]), &(energy_errors[0]), &(all_electron_fluxE3_errors[0]));

	gr_all_electron_flux.SetName("gr_all_electron_flux");
	gr_all_electron_flux_E3.SetName("gr_all_electron_flux_E3");
	gr_all_electron_flux.SetTitle("All Electron Flux");
	gr_all_electron_flux_E3.SetTitle("All Electron Flux");

	// Write to file
	if (verbose)
		std::cout << "\n\nFinalizing...";
	outFile.cd();
	auto fluxdir = outFile.mkdir("eFlux");
	fluxdir->cd();

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
#endif
}