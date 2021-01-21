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

extern const std::string uniqueOutFile(AnyOption &opt);

#endif