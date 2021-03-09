#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <vector>

#include "anyoption.h"

#pragma once

struct in_args
{
	std::string train_signal_input_list;
	std::string train_background_input_list;
	std::string test_signal_input_list;
	std::string test_background_input_list;
	std::string output_path;
	std::string config_dir;
	std::vector<std::string> learning_method;
	bool verbose = false;

	bool test_input_lists()
	{
		bool status = true;
		if (train_signal_input_list.empty() ||
			train_background_input_list.empty() ||
			test_signal_input_list.empty() ||
			test_background_input_list.empty() ||
			config_dir.empty())
			status = false;
		return status;
	}
	
	void ParseLearningMethods(std::string methods)
	{
		char sep = ':';
		size_t prev = 0;
		for(size_t ctr =0; ctr < methods.size(); ++ctr)
			if (methods[ctr]==sep)
			{
				learning_method.push_back(methods.substr(prev, ctr));
				prev = ctr+1;
			}
		if (!learning_method.size())
			learning_method.push_back(methods);
		else
			learning_method.push_back(methods.substr(prev, methods.size()-1));
		
	}	
};

extern void PrintArgs(in_args input_args);

#endif