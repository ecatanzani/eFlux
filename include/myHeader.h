#ifndef MYHEADER_H
#define MYHEADER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "anyoption.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TChain.h"

// DAMPESW includes

#include "DmpChain.h"


#pragma once

#define nGFitLoops 1000000
#define chiSQLLimit 0.000001

// DAMPE struct

#define DAMPE_bgo_nLayers 14


extern void eCore(
                    const std::string inputPath,
                    const std::string outputPath,
                    const bool verbose,
                    const bool pedantic,
                    const unsigned int lvTime,
                    const bool myAcceptance,
                    const std::string accInputPath,
                    AnyOption &opt
                );

extern void readInputTree(
                            const std::string inputPath,
                            std::vector<float> &dataValues,
                            TTree &dTree
                        );

extern void branchTree(TTree &myDataTree,std::vector<float> &dataValues);
extern const char* uniqueOutFile(const std::string outputPath,AnyOption &opt);
extern std::vector<float> createLogBinning(const double minValue,const double maxValue,const int nBins);
extern std::vector<float> createLinBinning(const double minValue,const double maxValue,const int nBins);

extern void buildAcceptance(
                                const std::string accInputPath,
                                const bool verbose,
                                const std::vector<float> &logEBins
                            );

extern DmpChain* aggregateEventsDmpChain(const std::string accInputPath,const bool verbose);
extern TChain* aggregateEventsTChain(const std::string accInputPath,const bool verbose);

extern std::string getListPath(const std::string accInputPath,const bool MC=false);

extern void buildFlux(
                        const std::string inputPath,
                        const unsigned int lvTime,
                        TFile &outFile,
                        const bool verbose,
                        const bool pedantic,
                        const std::string accInputPath,
                        const bool myAcceptance
                    );

extern bool chechFlags(
                        AnyOption &opt,
                        const std::string inputPath,
                        const std::string outputPath,
                        const unsigned int lvTime
                    );

extern void buildXtrlFlux(
                            std::vector<float> &eBins,
                            std::vector<float> &cBins,
                            const std::string inputPath,
                            const unsigned int lvTime,
                            TFile &outFile,
                            const bool verbose
                        );

extern void evLoop(
                    TH1D &inputHisto,
                    const std::string inputPath,
                    TFile &outFile,
                    const bool verbose,
                    const bool eClassifier=false,
                    const double xtrlCut=8.5
                );

extern void readAcceptance(TH1D &acceptance, TFile &outFile,const bool verbose);

extern TF1 fitGFactor(TH1D *gFactor,TFile &outFile, const bool verbose);

extern void tmpFit(
                    TDirectory* geoFactorDir,
                    TF1 &myFitter,
                    TH1D* gFactor,
                    const bool verbose,
                    const bool baseFit = true,
                    const unsigned int fitNumber = 0
                );

extern void setStartingParameters(
                                    const TF1 &oldFitter,
                                    TF1 &newFitter,
                                    const unsigned int nOldPars,
                                    const unsigned int nNewPars
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