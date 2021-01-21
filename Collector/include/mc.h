#ifndef MC_H
#define MC_H

#include <string>

#include "TFile.h"

#include "main.h"

#include "anyoption.h"

extern void mcCore(in_pars input_pars);

extern void mcLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool _VERBOSE,
	const std::string wd);

#endif