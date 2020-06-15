#include "acceptance_fit.h"
#include "TString.h"
#include "TGraph.h"
#include <memory>
#include <algorithm>
#include <numeric>

void tmpFit(
    TDirectory *accDir,
    TF1 &myFitter,
    TGraphAsymmErrors *acceptance,
    const bool verbose,
    const bool baseFit,
    const unsigned int fitNumber)
{
    // Set number of points to the TF1
    myFitter.SetNpx(100000);

    // Creating a TDirectory for the outputs of the tmo fit procedure on the Geometrical Factor

    TString tDirName = "logisticFunction";
    if (!baseFit)
    {
        tDirName += "_";
        tDirName += fitNumber;
    }
    auto logisticDir = accDir->mkdir(tDirName.Data());
    logisticDir->cd();

    // Getting the number of parameters from the TF1
    auto nPars = myFitter.GetNpar();

    // Creating Chi2,Prob and NDF histos
    TH1D chi2("chi2", "#chi^{2} - Geometric Factor Fit", 100, 0, 10);
    TH1D prob("prob", "Probability - Geometric Factor Fit", 100, 0, 10);
    TH1D ndf("ndf", "NDF - Geometric Factor Fit", 100, 0, 10);

    // vectors containing the parameters values
    std::vector<double> iters;
    std::vector<double> chi2iters;
    std::vector<std::vector<double>> fitPars(nPars);
    
    double chisq;   // Store chisq result

    // Start the fit procedure
    for (int fIdx = 0; fIdx < nGFitLoops; ++fIdx)
    {
        if (!fIdx)
            acceptance->Fit(myFitter.GetName(), "IR");    
        else
        {
            // Set starting parameter according to the previous fit
            for(auto pIdx = 0; pIdx < nPars; ++pIdx)
                myFitter.SetParameter(pIdx,myFitter.GetParameter(pIdx));
            acceptance->Fit(myFitter.GetName(), "IR");
            if (fabs(chisq - myFitter.GetChisquare()) < chiSQLLimit)
                break;
        }
        // Add the iteration number to iters vector
        iters.push_back(fIdx);
        for (int pIdx = 0; pIdx < nPars; ++pIdx)
            fitPars[pIdx].push_back(myFitter.GetParameter(pIdx));
        chi2.Fill(myFitter.GetChisquare());
        prob.Fill(myFitter.GetProb());
        ndf.Fill(myFitter.GetNDF());
        chi2iters.push_back(myFitter.GetChisquare());
        chisq = myFitter.GetChisquare();

        if (verbose)
        {
            std::cout << "\n\nIteration " << fIdx + 1;
            for (int pIdx = 0; pIdx < nPars; ++pIdx)
                std::cout << "\nPar " << pIdx << "\t-> " << myFitter.GetParameter(pIdx);
            std::cout << "\n Chi2 \t-> " << myFitter.GetChisquare();
            std::cout << "\n Prob \t-> " << myFitter.GetProb();
            std::cout << "\n NDF \t-> " << myFitter.GetNDF();
            std::cout << "\n Chi2/NDF \t-> " << myFitter.GetChisquare() / myFitter.GetNDF() << std::endl;
        }
    }

    std::vector<std::shared_ptr<TGraph>> grPar(nPars, nullptr);
    for (int pIdx = 0; pIdx < nPars; ++pIdx)
    {
        TString grName = "grPar_";
        TString grTitle = "grPar ";
        grName += pIdx;
        grTitle += pIdx;
        grPar[pIdx] = std::make_shared<TGraph>(iters.size(), &(iters[0]), &(fitPars[pIdx][0]));
        grPar[pIdx]->SetName(grName);
        grPar[pIdx]->SetTitle(grTitle);
        grPar[pIdx]->SetMarkerStyle(8);
    }

    // Create chi2 TGraph
    TGraph chi2gr(iters.size(), &(iters[0]), &(chi2iters[0]));
    chi2gr.SetName("chi2gr");
    chi2gr.SetTitle("#chi^{2} evolution");
    chi2gr.SetMarkerStyle(8);

    // Create pool TGraph
    std::vector<double> gr_points(acceptance->GetN());
    std::iota(std::begin(gr_points), std::end(gr_points), 0);
    std::vector<double> pull(acceptance->GetN());
    std::vector<double> fit_distance(acceptance->GetN());
    gr_point tmp_point;

    for(auto pIdx=0; pIdx<acceptance->GetN(); ++pIdx)
    {
        acceptance->GetPoint(pIdx, tmp_point.X, tmp_point.Y);
        pull[pIdx] = (tmp_point.Y-myFitter.Eval(tmp_point.X))/acceptance->GetErrorY(pIdx);
        fit_distance[pIdx] = tmp_point.Y-myFitter.Eval(tmp_point.X);
    }

    TGraph gr_pull(gr_points.size(), &(gr_points[0]), &(pull[0]));
    TGraph gr_distance(gr_points.size(), &(gr_points[0]), &(fit_distance[0]));
    gr_pull.SetName("gr_pull");
    gr_pull.SetTitle("pull");
    gr_distance.SetName("gr_distance");
    gr_distance.SetTitle("fit point distance");

    gr_pull.Write();
    gr_distance.Write();
    // Writing fitted Acceptance Factor
    acceptance->Write();

    // Writing chi2,prob and NDF histos
    chi2.Write();
    prob.Write();
    ndf.Write();
    chi2gr.Write();

    //Writing TGraphs
    for (int pIdx = 0; pIdx < nPars; ++pIdx)
        grPar[pIdx]->Write();

    // Returning to main Geometrical Factor directory
    accDir->cd();
}

inline void setStartingParameters(
    const TF1 &oldFitter,
    TF1 &newFitter,
    const unsigned int nOldPars,
    const unsigned int nNewPars)
{
    for (unsigned int pIdx = 0; pIdx < nOldPars; ++pIdx)
        newFitter.SetParameter(pIdx, oldFitter.GetParameter(pIdx));
}

TF1 fitAcceptance(TGraphAsymmErrors *finalAcceptance, TFile &outFile, const bool verbose)
{

    // Creating a TDirectory for the outputs of the fit procedure on the Geometrical Factor
    TDirectory *accDir = outFile.mkdir("Acceptance");
    accDir->cd();
    
    // Define start fitting TF1 - base logistic function
    TF1 fitter_0("fitter_0", logisticFunction_0, 1, 1e+4, 2);
    fitter_0.SetLineColor(kRed);
    tmpFit(accDir, fitter_0, finalAcceptance, verbose);
    
    // Define more accurate TF1s
    TF1 fitter_1("fitter_1", logisticFunction_1, 1, 1e+4, 4);
    setStartingParameters(fitter_0, fitter_1, fitter_0.GetNpar(), fitter_1.GetNpar());
    fitter_1.SetLineColor(kGreen + 1);
    tmpFit(accDir, fitter_1, finalAcceptance, verbose, false, 1);

    TF1 fitter_2("fitter_2", logisticFunction_2, 1, 1e+4, 6);
    setStartingParameters(fitter_1, fitter_2, fitter_1.GetNpar(), fitter_2.GetNpar());
    fitter_2.SetLineColor(kMagenta + 1);
    tmpFit(accDir, fitter_2, finalAcceptance, verbose, false, 2);
    
    TF1 fitter_3("fitter_3", logisticFunction_3, 1, 1e+4, 8);
    setStartingParameters(fitter_2, fitter_3, fitter_2.GetNpar(), fitter_3.GetNpar());
    fitter_3.SetLineColor(kCyan + 3);
    tmpFit(accDir, fitter_3, finalAcceptance, verbose, false, 3);

    /*
    TF1 fitter_4("fitter_4", logisticFunction_4, 1, 1e+4, 10);
    setStartingParameters(fitter_3, fitter_4, fitter_3.GetNpar(), fitter_4.GetNpar());
    fitter_4.SetLineColor(kBlue + 3);
    tmpFit(accDir, fitter_4, finalAcceptance, verbose, false, 4);
    
    TF1 fitter_5("fitter_5", logisticFunction_5, 1, 1e+4, 11);
    setStartingParameters(fitter_4, fitter_5, fitter_4.GetNpar(), fitter_5.GetNpar());
    fitter_5.SetLineColor(kOrange - 3);
    tmpFit(accDir, fitter_5, finalAcceptance, verbose, false, 5);
    
    TF1 fitter_6("fitter_6", logisticFunction_6, 1, 1e+4, 13);
    setStartingParameters(fitter_5, fitter_6, fitter_5.GetNpar(), fitter_6.GetNpar());
    fitter_6.SetLineColor(kCyan);
    tmpFit(accDir,fitter_6, finalAcceptance, verbose, false, 6);
    */
    /*
    TF1 fitter_7("fitter_7", logisticFunction_7, 1, 1e+4, 15);
    setStartingParameters(fitter_6, fitter_7, fitter_6.GetNpar(), fitter_7.GetNpar());
    fitter_6.SetLineColor(kGreen + 4);
    tmpFit(accDir, fitter_7, finalAcceptance, verbose, false, 7);
    */
   
    // Returning to main TFile directory
    outFile.cd();

    return fitter_3;
}