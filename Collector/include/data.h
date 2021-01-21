#ifndef DATA_H
#define DATA_H

#include <string>

#include "main.h"
#include "energy.h"
#include "anyoption.h"
#include "DmpNudContainer.h"
#include "DmpBgoContainer.h"
#include "DmpSimuTrajectory.h"
#include "DmpFilterContainer.h"

#include "TFile.h"
#include "TVector3.h"

extern void dataCore(in_pars input_pars);

extern void rawDataLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool _VERBOSE,
	const std::string wd);

#endif