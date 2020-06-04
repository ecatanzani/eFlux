#include "myHeader.h"
#include "acceptance_fit.h"

#include "TString.h"

template <typename InputDataType>
void tmpFit(
    TDirectory *accDir,
    TF1 &myFitter,
    InputDataType *acceptance,
    const bool verbose,
    const bool baseFit,
    const unsigned int fitNumber)
{
    // Set number of points to the TF1
    myFitter.SetNpx(1000000);

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
    TH1D chi2("chi2", "#chi^{2} - Geometric Factor Fit", 1000, 0, 10);
    TH1D prob("prob", "Probability - Geometric Factor Fit", 1000, 0, 1);
    TH1D ndf("ndf", "NDF - Geometric Factor Fit", 1000, 0, 40);

    // vectors containing the parameters values

#if 0
        std::vector<double> iters (nGFitLoops,.0);
        std::vector<double> chi2iters (nGFitLoops,.0);
        std::vector < std::vector <double> > fitPars (nPars, std::vector <double> (nGFitLoops));
#else
    std::vector<double> iters;
    std::vector<double> chi2iters;
    std::vector<std::vector<double>> fitPars(nPars);
#endif

    // Set default value for chisq
    double chisq;
    (*((long long *)&chisq)) = ~(1LL << 52);

    // Start the fit procedure
    for (int fIdx = 0; fIdx < nGFitLoops; ++fIdx)
    {
        acceptance->Fit(myFitter.GetName(), "IR");

        if (fabs(chisq - myFitter.GetChisquare()) < chiSQLLimit)
            break;

// Add the iteration number to iters vector
#if 0
            iters[fIdx] = fIdx;
#else
        iters.push_back(fIdx);
#endif

        for (int pIdx = 0; pIdx < nPars; ++pIdx)
#if 0
                fitPars[pIdx][fIdx] = myFitter.GetParameter(pIdx);
#else
            fitPars[pIdx].push_back(myFitter.GetParameter(pIdx));
#endif

            chi2.Fill(myFitter.GetChisquare());
        prob.Fill(myFitter.GetProb());
        ndf.Fill(myFitter.GetNDF());

#if 0
            chi2iters[fIdx] = myFitter.GetChisquare();
#else
        chi2iters.push_back(myFitter.GetChisquare());
#endif

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

    std::vector<TGraph *> grPar(nPars, nullptr);
    for (int pIdx = 0; pIdx < nPars; ++pIdx)
    {
        TString grName = "grPar_";
        TString grTitle = "grPar ";
        grName += pIdx;
        grTitle += pIdx;
        grPar[pIdx] = new TGraph(iters.size(), &(iters[0]), &(fitPars[pIdx][0]));
        grPar[pIdx]->SetName(grName);
        grPar[pIdx]->SetTitle(grTitle);
        grPar[pIdx]->SetMarkerStyle(8);
    }

    // Create chi2 TGraph
    TGraph chi2gr(iters.size(), &(iters[0]), &(chi2iters[0]));
    chi2gr.SetName("chi2gr");
    chi2gr.SetTitle("#chi^{2} evolution");
    chi2gr.SetMarkerStyle(8);

    // Writing fitted Geometrical Factor
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

    // Cleaning memory
    for (int pIdx = 0; pIdx < nPars; ++pIdx)
        grPar[pIdx]->Delete();
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

template <typename InputDataType>
TF1 fitAcceptance(InputDataType *finalAcceptance, TFile &outFile, const bool verbose)
{

    // Creating a TDirectory for the outputs of the fit procedure on the Geometrical Factor
    TDirectory *accDir = outFile.mkdir("Acceptance");
    accDir->cd();

    // Define start fitting TF1 - base logistic function
    TF1 fitter_0("fitter_0", logisticFunction_0, 1e-1, 1e+4, 4);
    fitter_0.SetLineColor(kRed);
    tmpFit<InputDataType>(accDir, fitter_0, finalAcceptance, verbose);

    // Define more accurate TF1s
    TF1 fitter_1("fitter_1", logisticFunction_1, 1e-1, 1e+4, 6);
    setStartingParameters(fitter_0, fitter_1, fitter_0.GetNpar(), fitter_1.GetNpar());
    fitter_1.SetLineColor(kGreen + 1);
    tmpFit<InputDataType>(accDir, fitter_1, finalAcceptance, verbose, false, 1);

    TF1 fitter_2("fitter_2", logisticFunction_2, 1e-1, 1e+4, 7);
    setStartingParameters(fitter_1, fitter_2, fitter_1.GetNpar(), fitter_2.GetNpar());
    fitter_2.SetLineColor(kMagenta + 1);
    tmpFit<InputDataType>(accDir, fitter_2, finalAcceptance, verbose, false, 2);

    TF1 fitter_3("fitter_3", logisticFunction_3, 1e-1, 1e+4, 9);
    setStartingParameters(fitter_2, fitter_3, fitter_2.GetNpar(), fitter_3.GetNpar());
    fitter_3.SetLineColor(kCyan + 3);
    tmpFit<InputDataType>(accDir, fitter_3, finalAcceptance, verbose, false, 3);

    TF1 fitter_4("fitter_4", logisticFunction_4, 1e-1, 1e+4, 11);
    setStartingParameters(fitter_3, fitter_4, fitter_3.GetNpar(), fitter_4.GetNpar());
    fitter_4.SetLineColor(kBlue + 3);
    tmpFit<InputDataType>(accDir, fitter_4, finalAcceptance, verbose, false, 4);

    TF1 fitter_5("fitter_5", logisticFunction_5, 1e-1, 1e+4, 13);
    setStartingParameters(fitter_4, fitter_5, fitter_4.GetNpar(), fitter_5.GetNpar());
    fitter_5.SetLineColor(kOrange - 3);
    tmpFit<InputDataType>(accDir, fitter_5, finalAcceptance, verbose, false, 5);

    TF1 fitter_6("fitter_6", logisticFunction_6, 1e-1, 1e+4, 15);
    setStartingParameters(fitter_5, fitter_6, fitter_5.GetNpar(), fitter_6.GetNpar());
    fitter_6.SetLineColor(kYellow + 2);
    tmpFit<InputDataType>(accDir, fitter_6, finalAcceptance, verbose, false, 6);

    TF1 fitter_7("fitter_7", logisticFunction_7, 1e-1, 1e+4, 18);
    setStartingParameters(fitter_6, fitter_7, fitter_6.GetNpar(), fitter_7.GetNpar());
    fitter_7.SetLineColor(kGreen + 4);
    tmpFit<InputDataType>(accDir, fitter_7, finalAcceptance, verbose, false, 7);

    /*
    TF1 fitter_8("fitter_8",logisticFunction_8,1e-1,1e+4,21);
    setStartingParameters(fitter_7,fitter_8,fitter_7.GetNpar(),fitter_8.GetNpar());
    fitter_8.SetLineColor(kCyan);
    tmpFit(accDir,fitter_8,acceptance,verbose,false,8);
    */

    // Returning to main TFile directory
    outFile.cd();

    return fitter_7;
}

template TF1 fitAcceptance(TH1F *finalAcceptance,TFile &outFile,const bool verbose);
template TF1 fitAcceptance(TH1D *finalAcceptance,TFile &outFile,const bool verbose);
template TF1 fitAcceptance(TGraphAsymmErrors *finalAcceptance,TFile &outFile,const bool verbose);
template void tmpFit(
    TDirectory *accDir,
    TF1 &myFitter,
    TH1D *acceptance,
    const bool verbose,
    const bool baseFit = true,
    const unsigned int fitNumber = 0);
template void tmpFit(
    TDirectory *accDir,
    TF1 &myFitter,
    TH1F *acceptance,
    const bool verbose,
    const bool baseFit = true,
    const unsigned int fitNumber = 0);
template void tmpFit(
    TDirectory *accDir,
    TF1 &myFitter,
    TGraphAsymmErrors *acceptance,
    const bool verbose,
    const bool baseFit = true,
    const unsigned int fitNumber = 0);