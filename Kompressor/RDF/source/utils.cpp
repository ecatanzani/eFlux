#include "utils.h"

void UpdateProcessStatus(
	const int evIdx,
	int &kStep,
	const int nevents)
{
	auto percentage = ((evIdx + 1) / (double)nevents) * 100;
	if (floor(percentage) != 0 && ((int)floor(percentage) % kStep) == 0)
	{
		std::cout << "\n"
				  << (int)percentage << " %\t | \tProcessed " << evIdx + 1 << " events / " << nevents;
		if ((int)percentage == 100)
			std::cout << "\n\n";
		kStep += 10;
	}
}