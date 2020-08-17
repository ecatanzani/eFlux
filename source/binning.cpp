#include "binning.h"

std::vector<float> createLogBinning(
	const double eMin,
	const double eMax,
	const std::size_t n_bins)
{
	std::vector<float> binning(n_bins + 1, 0);
	double log_interval = (log10(eMax) - log10(eMin)) / n_bins;
	for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
		binning[bIdx] = pow(10, log10(eMin) + bIdx * log_interval);

	return binning;
}

std::vector<float> createLinearBinning(
	float a,
	float b,
	std::size_t N)
{
	float h = (b - a) / static_cast<float>(N);
	std::vector<float> xs(N + 1);
	std::vector<float>::iterator x;
	float val;
	for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
		*x = val;

	return xs;
}