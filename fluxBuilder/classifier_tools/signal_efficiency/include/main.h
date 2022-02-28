#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <memory>

#include "anyoption.h"

#pragma once

struct in_args
{
    const char* energy_config_file;
	std::string input_list;
	std::string output_path;
	std::string eff_corr_function;
	bool verbose {false};
	unsigned int threads {1};
};

#endif