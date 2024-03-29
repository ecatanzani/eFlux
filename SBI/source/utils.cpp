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
		kStep += 10;
	}
}

const std::string uniqueOutFile(AnyOption &opt)
{
	std::time_t ctime = std::time(0);
	std::stringstream fPath;
	if (opt.getValue("outputDir") || opt.getValue('d'))
		fPath << opt.getValue('d') << "/SBIOutFile_" << ctime << ".root";
	else if (opt.getValue("output") || opt.getValue('o'))
		fPath << opt.getValue('o');
	else
		fPath << "SBIOutFile_" << ctime << ".root";

	return fPath.str();
}