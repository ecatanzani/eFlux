#ifndef SHOWGAUS_H
#define SHOWGAUS_H

#include "main.h"

#include <string>
#include <vector>
#include <memory>
#include <iostream>

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <ROOT/RDataFrame.hxx>

extern void showgaus(in_args input_args);
extern void resumegaus(in_args input_args);

#endif