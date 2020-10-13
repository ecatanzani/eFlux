#ifndef MYHEADER_H
#define MYHEADER_H

#include <iostream>
#include <string>
#include <vector>

#include "anyoption.h"

#pragma once

extern void eCore(
    const std::string inputPath,
    const std::string outputPath,
    const bool verbose,
    const bool pedantic,
    const bool rawdata_flag,
	const bool skimmed_flag,
	const bool ntuple_flag,
    AnyOption &opt,
    const std::string wd);

extern const std::string uniqueOutFile(
    const std::string outputPath,
    AnyOption &opt);

#if 0
extern const std::string uniqueTupleOutFile(
    AnyOption &opt, 
    const int year, 
    const int month,
    const int emin,
    const int emax);
#else
extern const std::string uniqueTupleOutFile(
    AnyOption &opt, 
    const std::string year, 
    const std::string month);
#endif

extern bool chechFlags(
    AnyOption &opt,
    const std::string inputPath,
    const std::string outputPath,
    const unsigned int lvTime);

extern void generateFinalGraph(
    const bool verbose,
    const bool pedantic,
    const std::string outputPath,
    const std::string complete_histo_path,
    const std::string wd);

extern void generateDataFinalGraph(
    const bool verbose,
    const bool pedantic,
    const std::string outputPath,
    const std::string complete_histo_path,
    const std::string wd);

#endif