#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>

#include "anyoption.h"

#ifdef _DEBUG 
	#define _RooFitVerbosity true
#else
	#define _RooFitVerbosity false
#endif

#define _RooFitClamping true
struct deps_paths
{
	std::string data_path;
	std::string electron_path;
	std::string proton_path;
};

#endif