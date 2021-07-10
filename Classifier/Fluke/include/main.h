#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <vector>

#include "anyoption.h"

#pragma once

struct in_args
{
	std::string input_list;
	std::string output_path;
	std::string reg_fit_path;
	std::string lambda_path;
	std::string wd;

	bool verbose = false;
	bool mc_flag = false;
	bool split;
	bool split_test;
	bool signal;

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
		split = true;
		if (!set_name.empty())
		{
			if (!strcmp(set_name.c_str(), "t"))
				split_test = false;
			else if (!strcmp(set_name.c_str(), "T"))
				split_test = true;
			else
			{
				std::cerr << "\n\nError parsing tmva no-split option... exiting\n\n";
				exit(100);
			}
		}
	}
};

#endif