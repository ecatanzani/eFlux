#include "myHeader.h"

#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "TString.h"

double logisticFunction(double *x, double *par)
{
    double xx =x[0];
    double logistic = par[0] + ( par[1] ) / ( 1 + TMath::Exp( - par[2] * ( par[3]*TMath::Power(xx,2) + par[4]*xx +par[5] ) ) );
    return logistic;
}

void fitGFactor(TH1D *gFactor, TFile &outFile, const bool verbose)
{
    int nPars = 6;
    
    //TF1 gFitter("gFitter",logisticFunction,0,1e+4,nPars);
    TF1 gFitter("gFitter","[0] + ( [1] ) / ( 1 + TMath::Exp( - [2] * ( [3]*TMath::Power(x,2) + [4]*x +[5] ) ) ) ",0,1e+4);

    // vectors containing the parameters values
    std::vector<double> iters(nGFitLoops,.0); 
    std::vector < std::vector <double> > fitPars(nPars, std::vector <double> (nGFitLoops));

    for(int fIdx=0; fIdx<nGFitLoops; ++fIdx)
    {
        // Add the iteration number to iters vector
        iters[fIdx] = fIdx;
        
        gFactor->Fit("gFitter","I");

        for(int pIdx=0; pIdx<nPars; ++pIdx)
            fitPars[pIdx][fIdx] = gFitter.GetParameter(pIdx);
        
        if(verbose)
        {
            std::cout << "\n\nIteration " << fIdx+1;
            for(int pIdx=0; pIdx<nPars; ++pIdx)
                std::cout << "\nPar " << pIdx <<"\t-> " << gFitter.GetParameter(pIdx);
        }
    }
    
    TGraph* grPar[nPars];
    for(int pIdx=0; pIdx<nPars; ++pIdx)
    {
        TString grName = "grPar_";
        TString grTitle = "grPar ";
        grName += pIdx;
        grTitle += pIdx;
        grPar[pIdx] = new TGraph(nGFitLoops, & (iters[0]), & (fitPars[pIdx][0]));
        grPar[pIdx]->SetName(grName);
        grPar[pIdx]->SetTitle(grTitle);
    }
    
    // Creating a TDirectory for the outputs of the fit procedure on the Geometrical Factor
    TDirectory *geoFactorDir = outFile.mkdir("geoFactor");
    geoFactorDir->cd();

    // Writing fitted Geometrical Factor
    gFactor->Write();

    //Writing TGraphs
    for(int pIdx=0; pIdx<nPars; ++pIdx)
        grPar[pIdx]->Write();

    // Returning to main dir on the output TFile
    outFile.cd();

    // Cleaning memory
    for(int pIdx=0; pIdx<nPars; ++pIdx)
        grPar[pIdx]->Delete();
}