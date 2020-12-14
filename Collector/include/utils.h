#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <math.h>
#include <ctime>
#include <memory>

#include "anyoption.h"

extern void UpdateProcessStatus(
	const int evIdx,
	int &kStep,
	const int nevents);

extern const std::string uniqueOutFile(
	const std::string outputPath,
	AnyOption &opt);

extern std::shared_ptr<ofstream> ev_logger(
	const std::string outputPath,
	AnyOption &opt);

#endif