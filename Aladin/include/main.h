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
	std::string best_lambda_tree_path;
	bool verbose = false;

	// Tasks
	bool mc_flag = false;
	bool gaussianize = false;
	bool fit = false;
	bool loglikelihood = false;
	unsigned int energybin = 999;
	unsigned int threads = 1;

	const bool check_input()
	{	
		bool status = false;
		if (loglikelihood || fit)
		{
			if (!wd.empty() && !input_list.empty() && !output_path.empty() && energybin!=999) 
				status = true;
		}
		else if (gaussianize)
		{
			if (!wd.empty() && !input_list.empty() && !output_path.empty() && !regularize_tree_path.empty()) 
				status = true;
		}
		else
		{
			if (!best_lambda_tree_path.empty() && !regularize_tree_path.empty())
				status = true;
		}

		return status;
	}
};

#endif