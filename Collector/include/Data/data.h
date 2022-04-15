#ifndef DATA_H
#define DATA_H

#include <string>

#include "main.h"
#include "energy.h"
#include "anyoption.h"
#include "DmpSimuTrajectory.h"
#include "Dmp/DmpNudContainer.h"
#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpFilterContainer.h"

#include "TFile.h"
#include "TVector3.h"

extern void dataCore(in_pars input_pars);

extern void rawDataLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool verbose,
	const std::string wd,
	const std::string stk_corrections);

#endif