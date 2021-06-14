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
	std::string fit_tree_path;
	bool _VERBOSE = false;

	// Tasks
	bool mc_flag = false;
	bool tmva_set = false;
	bool signal;
	bool no_split = false;
	bool no_split_test;
	bool gaussianize = false;
	bool study_gaussianized = false;
	bool fit_gaussianized = false;
	bool loglikelihood = false;
	unsigned int likelihood_energybin = 0;
	unsigned int threads = 1;

	void SetSeType(std::string set_name)
	{
		if (!set_name.empty())
		{
			if (!strcmp(set_name.c_str(), "s") || !strcmp(set_name.c_str(), "S"))
				signal = true;
			else if (!strcmp(set_name.c_str(), "b") || !strcmp(set_name.c_str(), "B"))
				signal = false;
			else
			{
				std::cerr << "\n\nError parsing tmva-set option... exiting\n\n";
				exit(100);
			}
		}
	}
	void SetNSeType(std::string set_name)
	{
		no_split = true;
		if (!set_name.empty())
		{
			if (!strcmp(set_name.c_str(), "t"))
				no_split_test = false;
			else if (!strcmp(set_name.c_str(), "T"))
				no_split_test = true;
			else
			{
				std::cerr << "\n\nError parsing tmva no-split option... exiting\n\n";
				exit(100);
			}
		}
	}
};

#endif