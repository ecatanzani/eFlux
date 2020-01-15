#ifndef MYHEADER_H
#define MYHEADER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"

#pragma once

extern void eCore(
                    const std::string inputPath,
                    const std::string outputPath,
                    const bool verbose
                );

extern void readInputTree(const std::string inputPath,std::vector<double> &dataValues,TTree* dTree);
extern void branchTree(TTree &myDataTree,std::vector<double> &dataValues);
extern const char* uniqueTree(const std::string outputPath);
extern void uniqueTree(TTree *dTree,const TFile &outFile);
extern std::vector<double> createLogBinning(const double minValue,const double maxValue,const int nBins);
extern void buildFlux(TFile &outFile);
extern void buildXtrlFlux(TFile &outFile,std::vector<double> &logBins);
extern void evLoop(TH2D &inputHisto,const std::string inputPath,const bool eClassifier);

#endif