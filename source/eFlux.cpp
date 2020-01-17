#include "myHeader.h"

void buildFlux(const std::string inputPath,const unsigned int lvTime)
{
    
    double minValue = 3e-10;
    double maxValue = 5e+10;
    const int nBins = 10;

    std::vector<double> logEBins = createLogBinning(minValue,maxValue,nBins);
    
    minValue = 10;
    maxValue = 1e+5;
    
    std::vector<double> eCounts = createLinBinning(minValue,maxValue,nBins);

    buildXtrlFlux(logEBins,eCounts,inputPath,lvTime);

}

void buildXtrlFlux(
                    std::vector<double> &eBins,
                    std::vector<double> &cBins,
                    const std::string inputPath,
                    const unsigned int lvTime
                )
{
    /*
        All-Electron flux using xtrl as classifier
    */

    TH1D eCounts("eCounts","All Electron counts - xtrl classifier",eBins.size(),&(eBins[0]));
    TH1D acceptance;
    TH1D* eFlux = nullptr;

    evLoop(eCounts,inputPath,true);
    buildAcceptance(acceptance);

    eFlux = (TH1D*)eCounts.Clone("eFlux");
    eFlux->Sumw2();
    eFlux->Reset();

    //Building flux
    for(int bIdx=1; bIdx<=eFlux->GetXaxis()->GetNbins(); ++bIdx)
        eFlux->SetBinContent(bIdx,eCounts.GetBinContent(bIdx)/acceptance.GetBinContent(bIdx));
    
    //Scale flux for the live-time
    eFlux->Scale(1/(double)lvTime);

    eFlux->Write();

}