#include "myHeader.h"
#include "acceptance.h"
#include "flux.h"
#include "dataLoop.h"

void buildXtrlFlux(
    const std::vector<float> &eBins,
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose,
    const std::string accInputPath)
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
        accInputPath);

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