#ifndef MC_H
#define MC_H

#include <string>

#include "TFile.h"

#include "anyoption.h"

extern void mcCore(
	const std::string inputPath,
	const std::string outputPath,
	const bool _VERBOSE,
	const bool pedantic,
	AnyOption &opt,
	const std::string wd);

extern void mcLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool verbose,
	const std::string wd);

#endif