#include "myHeader.h"

#include "TDirectory.h"

void buildFlux(
                const std::string inputPath,
                const unsigned int lvTime,
                TFile &outFile,
                const bool verbose
            )
{
    
    double minValue = 3e-10;
    double maxValue = 5e+10;
    const int nBinsE = 10;
    const int nBinsC = 10;

    std::vector<float> logEBins = createLogBinning(minValue,maxValue,nBinsE);
    
    if(verbose)
    {
        std::cout << "\nEnergy log binning...";
        for(unsigned int idx=0; idx < logEBins.size(); ++idx)
            std::cout << "\n" << logEBins[idx];
    }

    minValue = 0;
    maxValue = 1e+5;
    
    std::vector<float> eCounts = createLinBinning(minValue,maxValue,nBinsC);
    
    if(verbose)
    {
        std::cout << "\nCounts linear binning...";
        for(unsigned int idx=0; idx < eCounts.size(); ++idx)
            std::cout << "\n" << eCounts[idx];
    }

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

    evLoop(eCounts,inputPath,outFile,true);
    buildAcceptance(acceptance,outFile,verbose);
    
    eFlux = (TH1D*)eCounts.Clone("eFlux");
    eFlux->Sumw2();
    eFlux->Reset();

    //Building flux
    for(int bIdx=1; bIdx<=eFlux->GetXaxis()->GetNbins(); ++bIdx)
        eFlux->SetBinContent(bIdx,eCounts.GetBinContent(bIdx)/acceptance.GetBinContent(bIdx));
    
    //Scale flux for the live-time
    eFlux->Scale(1/(double)lvTime);

    // Creating a TDirectory for the flux histo
    TDirectory *fDir = outFile.mkdir("Flux");
    fDir->cd();

    eFlux->Write();

    // Returning to main dir on the output TFile
    outFile.cd();
}