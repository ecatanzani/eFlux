#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include "anyoption.h"

#pragma once

struct in_args
{
    // Vars
	std::string input_list;
	std::string output_path;
	bool verbose = false;

	// Tasks
	bool mc_flag = false;
	unsigned int threads = 1;
};

#endif