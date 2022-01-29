#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <memory>

#include "anyoption.h"

#pragma once

struct in_args
{
    std::string input_list;
	const char* energy_config_file;
	const char* bdt_config_file;
	const char* acceptance_file;
	unsigned int energy_bin;
	double exposure;
	std::string output_path;
	std::string learning_method;
	bool verbose {false};
};

#endif