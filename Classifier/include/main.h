#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>

#include "anyoption.h"

#pragma once

struct in_args
{
    std::string train_signal_path;
	std::string train_background_path;
	std::string test_signal_path;
	std::string test_background_path;
	std::string output_path;
	std::string learning;
	bool verbose = false;
	bool test_with_data = false;
};

extern void PrintArgs(in_args input_args);

#endif