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
};

#endif