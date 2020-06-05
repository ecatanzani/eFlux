#ifndef ACCEPTANCE_FIT_H
#define ACCEPTANCE_FIT_H

#include "TF1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TDirectory.h"

#define nGFitLoops 1000000
#define chiSQLLimit 0.000001

template <typename InputDataType> extern TF1 fitAcceptance(
    InputDataType *gFactor, 
    TFile &outFile, 
    const bool verbose);

extern template TF1 fitAcceptance(TH1F *finalAcceptance,TFile &outFile,const bool verbose);
extern template TF1 fitAcceptance(TH1D *finalAcceptance,TFile &outFile,const bool verbose);
extern template TF1 fitAcceptance(TGraphAsymmErrors *finalAcceptance,TFile &outFile,const bool verbose);

template <typename InputDataType> extern void tmpFit( 
    TDirectory *accDir,
    TF1 &myFitter,
    InputDataType *gFactor,
    const bool verbose,
    const bool baseFit = true,
    const unsigned int fitNumber = 0);

extern template void tmpFit(
    TDirectory *accDir,
    TF1 &myFitter,
    TH1D *acceptance,
    const bool verbose,
    const bool baseFit = true,
    const unsigned int fitNumber = 0);

extern template void tmpFit(
    TDirectory *accDir,
    TF1 &myFitter,
    TH1F *acceptance,
    const bool verbose,
    const bool baseFit = true,
    const unsigned int fitNumber = 0);

extern template void tmpFit(
    TDirectory *accDir,
    TF1 &myFitter,
    TGraph *acceptance,
    const bool verbose,
    const bool baseFit = true,
    const unsigned int fitNumber = 0);

// Fitting functions
extern double logisticFunction_0(double *x, double *par);
extern double logisticFunction_1(double *x, double *par);
extern double logisticFunction_2(double *x, double *par);
extern double logisticFunction_3(double *x, double *par);
extern double logisticFunction_4(double *x, double *par);
extern double logisticFunction_5(double *x, double *par);
extern double logisticFunction_6(double *x, double *par);
extern double logisticFunction_7(double *x, double *par);
extern double logisticFunction_8(double *x, double *par);

#endif