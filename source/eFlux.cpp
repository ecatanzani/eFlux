#include "myHeader.h"
#include "acceptance.h"

#include "TDirectory.h"

void buildFlux(
                const std::string inputPath,
                const unsigned int lvTime,
                TFile &outFile,
                const bool verbose,
                const bool pedantic,
                const std::string accInputPath,
                const bool myAcceptance
            )
{
    
    // max and min value of the energy (in GeV)
    double minValue = 1;
    double maxValue = 1e+4;
    const int nBinsE = 1000;
    const int nBinsC = 1000;

    std::vector<float> logEBins = createLogBinning(minValue,maxValue,nBinsE);

    if(pedantic)
    {
        std::cout << "\nEnergy log binning..." << std::scientific;
        for(unsigned int idx=0; idx < logEBins.size(); ++idx)
            std::cout << "\n" << logEBins[idx];
        std::cout << std::defaultfloat;
    }
    
    minValue = 0;
    maxValue = 1e+5;
    
    std::vector<float> eCounts = createLinBinning(minValue,maxValue,nBinsC);
    
    if(pedantic)
    {
        std::cout << "\nCounts linear binning...";
        for(unsigned int idx=0; idx < eCounts.size(); ++idx)
            std::cout << "\n" << eCounts[idx];
    }

    // Building acceptance
    if(myAcceptance)
        buildAcceptance(accInputPath,verbose,logEBins);

    buildXtrlFlux(logEBins,eCounts,inputPath,lvTime,outFile,verbose);

}

void buildXtrlFlux(
                    std::vector<float> &eBins,
                    std::vector<float> &cBins,
                    const std::string inputPath,
                    const unsigned int lvTime,
                    TFile &outFile,
                    const bool verbose
                )
{
    /*
        All-Electron flux using xtrl as classifier
    */

    TH1D eCounts("eCounts","All Electron counts - xtrl classifier",eBins.size()-1,&(eBins[0]));
    TH1D acceptance;
    TH1D* eFlux = nullptr;

    evLoop(eCounts,inputPath,outFile,verbose,true);
    readAcceptance(acceptance,outFile,verbose);
    
    eFlux = (TH1D*)eCounts.Clone("eFlux");
    eFlux->Sumw2();
    eFlux->Reset();

    //Building flux
    for(int bIdx=1; bIdx<=eFlux->GetXaxis()->GetNbins(); ++bIdx)
        if(acceptance.GetBinContent(bIdx))
            eFlux->SetBinContent(bIdx,eCounts.GetBinContent(bIdx)/acceptance.GetBinContent(bIdx));

    //Scale flux for the live-time
    eFlux->Scale(1/(double)lvTime);

    // Creating a TDirectory for the flux histo
    TDirectory *fDir = outFile.mkdir("Flux");
    fDir->cd();

    eCounts.Write();
    eFlux->Write();

    // Returning to main dir on the output TFile
    outFile.cd();
}