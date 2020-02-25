#ifndef ACCEPTANCEFIT_H
#define ACCEPTANCEFIT_H

#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TDirectory.h"

extern void setStartingParameters(
                                    const TF1 &oldFitter,
                                    TF1 &newFitter,
                                    const unsigned int nOldPars,
                                    const unsigned int nNewPars
                                );

extern TF1 fitGFactor(TH1D *gFactor,TFile &outFile, const bool verbose);

extern void tmpFit(
                    TDirectory* geoFactorDir,
                    TF1 &myFitter,
                    TH1D* gFactor,
                    const bool verbose,
                    const bool baseFit = true,
                    const unsigned int fitNumber = 0
                );

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