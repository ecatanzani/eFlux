#ifndef BINNING_H
#define BINNING_H

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

extern std::vector<float> createLogBinning(
	const double eMin, 
	const double eMax, 
	const size_t n_bins);

extern std::vector<float> createLinearBinning(
	float a, 
	float b, 
	std::size_t N);

#endif