#ifndef DATA_H
#define DATA_H

#include <string>
#include <fstream>
#include <memory>

#include "anyoption.h"

#include "TFile.h"

extern void dataCore(
	const std::string inputPath,
	const std::string outputPath,
	const bool verbose,
	const bool pedantic,
	const bool edfilter_flag,
	AnyOption &opt,
	const std::string wd);

extern void rawDataLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool _VERBOSE,
	std::shared_ptr<ofstream> evlogger,
	const std::string wd);

#endif