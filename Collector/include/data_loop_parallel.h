#ifndef DATA_LOOP_PARALLEL_H
#define DATA_LOOP_PARALLEL_H

#include <string>

#include <ROOT/TSeq.hxx>
#include "TFile.h"
#include "TROOT.h"

extern void skimmedDataLoop_parallel(
	const std::string inputPath,
	TFile &outFile,
	const bool verbose,
	const bool skimmed,
	const std::string wd);

#endif