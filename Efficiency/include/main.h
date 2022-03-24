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
	std::string input_list;
	std::string bdt_config_file;
	std::string output_path;
	std::string energy_config_file;
	std::string bdt_learning_method;
	std::string cosine_regularize_path;
	std::string box_cox_regularize_path;
	std::string eff_corr_function;
	bool verbose {false};

	// Tasks
	unsigned int threads = 1;
};

#endif