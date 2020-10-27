#ifndef MC_H
#define MC_H

#include <string>

#include "TFile.h"

#include "anyoption.h"

extern void mcCore(
	const std::string inputPath,
	const bool _VERBOSE,
	const bool pedantic,
	const std::string outputPath,
	AnyOption &opt,
	const std::string wd);

extern void mcLoop(
	const std::string inputPath,
	const bool verbose,
	TFile &outFile,
	const std::string wd);

#endif