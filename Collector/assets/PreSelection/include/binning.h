#ifndef BINNING_H
#define BINNING_H

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

extern std::vector<double> createLogBinning(
	const double eMin, 
	const double eMax, 
	const size_t n_bins);

extern std::vector<double> createLinearBinning(
	double a, 
	double b, 
	std::size_t N);

#endif