#ifndef MYHEADER_H
#define MYHEADER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"

#pragma once

extern void eCore(
                    const std::string inputPath,
                    const std::string outputPath,
                    const bool verbose
                );

extern void readInputTree(const std::string inputPath,std::vector<double> &dataValues,TTree* dTree);
extern void branchTree(TTree &myDataTree,std::vector<double> &dataValues);
extern const char* uniqueOutFile(const std::string outputPath);
extern std::vector<double> createLogBinning(const double minValue,const double maxValue,const int nBins);
extern std::vector<double> createLinBinning(const double minValue,const double maxValue,const int nBins);
extern void buildFlux(TFile &outFile,const std::string inputPath);
extern void buildXtrlFlux(TFile &outFile,std::vector<double> &eBins,std::vector<double> &cBins,const std::string inputPath);

extern void evLoop(
                    TFile &outFile,
                    TH1D &inputHisto,
                    const std::string inputPath,
                    const bool eClassifier=false,
                    const double xtrlCut=8.5
                );

#endif