#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>

#include "anyoption.h"

#pragma once

struct in_args
{
	std::string train_signal_input_list;
	std::string train_background_input_list;
	std::string test_signal_input_list;
	std::string test_background_input_list;
	std::string output_path;
	std::string learning_method;
	bool verbose = false;
	bool test_with_data = false;

	bool test_input_lists()
	{
		bool status = true;
		if (train_signal_input_list.empty() ||
			train_background_input_list.empty() ||
			test_signal_input_list.empty() ||
			test_background_input_list.empty())
			status = false;
		return status;
	}
};

extern void PrintArgs(in_args input_args);

#endif