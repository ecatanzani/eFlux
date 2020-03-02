#include "myHeader.h"
#include "acceptance.h"
#include "flux.h"
#include "dataLoop.h"

#include "TDirectory.h"

void buildFlux(
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose,
    const bool pedantic,
    const std::string accInputPath,
    const bool myAcceptance)
{
    // Reading flux config from file
    flux_conf fluxParams;
    load_flux_struct(fluxParams);

    auto logEBins = createLogBinning(fluxParams);

    if (pedantic)
    {
        std::cout << "\nEnergy log binning..." << std::scientific;
        for (auto it = logEBins.begin(); it != logEBins.end(); ++it)
            std::cout << "\n"
                      << *it;
        std::cout << std::defaultfloat;
    }

    // Building acceptance
    if (myAcceptance)
        buildAcceptance(
            accInputPath,
            verbose,
            logEBins,
            outFile);

    buildXtrlFlux(
        logEBins,
        inputPath,
        lvTime,
        outFile,
        verbose,
        myAcceptance);
}

/*
    All-Electron flux using xtrl as classifier
*/
void buildXtrlFlux(
    const std::vector<float> &eBins,
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose,
    const bool myAcceptance)
{
    auto e_dataCounts =
        evLoop(
            eBins,
            inputPath,
            outFile,
            verbose,
            true);

    auto acceptance_tf1 = readAcceptance(
        outFile,
        verbose,
        myAcceptance);

    /*
    TH1D *eFlux = nullptr;
    eFlux = (TH1D *)eCounts.Clone("eFlux");
    eFlux->Sumw2();
    eFlux->Reset();
    */
   /*
    //Building flux
    for (int bIdx = 1; bIdx <= eFlux->GetXaxis()->GetNbins(); ++bIdx)
        if (acceptance.GetBinContent(bIdx))
            eFlux->SetBinContent(bIdx, eCounts.GetBinContent(bIdx) / acceptance.GetBinContent(bIdx));

    //Scale flux for the live-time
    eFlux->Scale(1 / (double)lvTime);

    // Creating a TDirectory for the flux histo
    TDirectory *fDir = outFile.mkdir("Flux");
    fDir->cd();

    eCounts.Write();
    eFlux->Write();

    // Returning to main dir on the output TFile
    outFile.cd();
    */
}

std::vector<float> createLogBinning(const flux_conf flux_params)
{
    std::vector<float> binning(flux_params.n_energy_bins + 1, 0);

    const double minLog = std::log(flux_params.min_energy);
    const double maxLog = std::log(flux_params.max_energy);
    const double logIncrement = (maxLog - minLog) / flux_params.n_energy_bins;
    double value = flux_params.min_energy;
    double logValue = minLog;

    for (auto it = binning.begin(); it != (binning.end() - 1); ++it)
    {
        auto bIdx = std::distance(binning.begin(), it);
        binning[bIdx] = value;
        logValue += logIncrement;
        value = std::exp(logValue);
    }
    binning[flux_params.n_energy_bins] = value;

    return binning;
}