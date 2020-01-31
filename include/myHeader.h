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

#pragma once

#define nGFitLoops 100

extern void eCore(
                    const std::string inputPath,
                    const std::string outputPath,
                    const bool verbose,
                    const unsigned int lvTime,
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

extern void buildFlux(
                        const std::string inputPath,
                        const unsigned int lvTime,
                        TFile &outFile,
                        const bool verbose
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

extern void buildAcceptance(TH1D &acceptance, TFile &outFile,const bool verbose);

extern void fitGFactor(TH1D *gFactor,TFile &outFile, const bool verbose);
extern double logisticFunction(double *x, double *par);

#endif