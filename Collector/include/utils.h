#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <cmath>
#include <math.h>
#include <ctime>
#include <sstream>

#include "anyoption.h"

extern void UpdateProcessStatus(
	const int evIdx,
	int &kStep,
	const int nevents);

extern const std::string uniqueOutFile(
	const std::string outputPath,
	AnyOption &opt);

#endif