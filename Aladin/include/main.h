#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <memory>

#include "anyoption.h"

#pragma once

struct in_args
{
    // Vars
    std::string wd;
	std::string input_list;
	std::string output_path;
	std::string regularize_tree_path;
	bool verbose = false;

	// Tasks
	bool mc_flag = false;
	bool gaussianize = false;
	bool study_gaussianized = false;
	bool fit_gaussianized = false;
	bool loglikelihood = false;
	unsigned int likelihood_energybin = 0;
	unsigned int threads = 1;

	const bool check_paths()
	{	
		if (loglikelihood)
		{
			if (!wd.empty() && !input_list.empty() && !output_path.empty()) return true;
			else return false;
		}
		else
		{
			if (!wd.empty() && !input_list.empty() && !output_path.empty() && !regularize_tree_path.empty()) return true;
			else return false;
		}
	}
};

#endif