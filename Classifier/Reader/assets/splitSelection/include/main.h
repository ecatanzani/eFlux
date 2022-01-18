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
	std::string output_directory;
	const char* energy_config_file;
	bool verbose {false};
	unsigned int threads {1};
};

#endif