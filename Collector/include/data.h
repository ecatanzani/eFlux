#ifndef DATA_H
#define DATA_H

#include <string>

#include "anyoption.h"

#include "TFile.h"

extern void dataCore(
	const std::string inputPath,
	const std::string outputPath,
	const bool verbose,
	const bool pedantic,
	const bool rawdata_flag,
	const bool ntuple_flag,
	AnyOption &opt,
	const std::string wd);

extern void rawDataLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool _VERBOSE,
	const std::string wd);

extern void tupleLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool _VERBOSE,
	const std::string wd);

#endif