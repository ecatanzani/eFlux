#include "myHeader.h"

#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TDirectory.h"

double logisticFunction(double *x, double *par)
{
    float xx =x[0];
    double logistic = par[0]/( 1 + TMath::Exp(-par[1] * (xx - par[2])) );
    return logistic;
}

void fitGFactor(TH1D *gFactor, TFile &outFile)
{
    int nPars = 3;
    TF1 gFitter("gFitter",logisticFunction,0,1e+4,nPars);
    
    // vectors containing the parameters values
    std::vector<double> iters(nGFitLoops,.0); 
    std::vector < std::vector <double> > fitPars(nPars, std::vector <double> (nGFitLoops));

    for(int fIdx=0; fIdx<nGFitLoops; ++fIdx)
    {
        // Add the iteration number to iters vector
        iters[fIdx] = fIdx;
        
        for(int pIdx=0; pIdx<nPars; ++pIdx)
        {
            // Set initial values for the fit.
            // In case of first iteration manually set initial values, otherwise get the fitted values obtained at the previous iteration
            if(fIdx==0 && pIdx==0)
            {
                gFitter.SetParameter(0,0.3);
                gFitter.SetParameter(1,1);
                gFitter.SetParameter(2,4e+3);
            }
            else
                for(int pIdx=0; pIdx<nPars; ++pIdx)
                    gFitter.SetParameter(pIdx,gFitter.GetParameter(pIdx));
            
            // Fixing 1 parameter on 3 on the fit model...            
            for(int fIdx=0; fIdx<3; ++fIdx)
                if(fIdx!=pIdx)
                    gFitter.FixParameter(fIdx,gFitter.GetParameter(fIdx));
            
            // Fit the model to the MC
            gFactor->Fit("gFitter","VIWM0");

            // Release fixed parameter
            for(int fIdx=0; fIdx<3; ++fIdx)
                if(fIdx!=pIdx)
                    gFitter.ReleaseParameter(fIdx);
            
        }

        // Add the fit parameters value to fitPars vector
        for(int pIdx=0; pIdx<nPars; ++pIdx)
            fitPars[pIdx][fIdx] = gFitter.GetParameter(pIdx);
        
    }
    
    TGraph grPar0(nGFitLoops, & (iters[0]), & (fitPars[0][0]));
    TGraph grPar1(nGFitLoops, & (iters[0]), & (fitPars[1][0]));
    TGraph grPar2(nGFitLoops, & (iters[0]), & (fitPars[2][0]));
    
    // Creating a TDirectory for the outputs of the fit procedure on the Geometrical Factor
    TDirectory *geoFactorDir = outFile.mkdir("geoFactor");
    geoFactorDir->cd();

    // Writing fitted Geometrical Factor
    gFactor->Write();

    //Writing TGraphs
    grPar0.Write();
    grPar1.Write();
    grPar2.Write();

    // Returning to main dir on the output TFile
    outFile.cd();
}