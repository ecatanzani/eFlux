#ifndef SHOWGAUS_H
#define SHOWGAUS_H

#include "main.h"

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <ROOT/RDataFrame.hxx>

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>

#define bgolayers 14
#define nenergybin 50

extern std::string parse_input_file(const std::string input_list);
extern std::shared_ptr<TChain> getchain(
    const std::string filelist, 
    const bool mc = false,
    const bool verbose = true);
extern std::vector<std::vector<std::shared_ptr<TH1D>>> getrmslayerhistos();
extern void showgaus(in_args input_args);

#endif