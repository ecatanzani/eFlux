#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <memory>

#include "anyoption.h"

#pragma once

struct in_args {
    // Vars
	std::string input_list;
	std::string output_path;
	std::string working_dir;
	bool split {true};
	unsigned int threads {1};
	bool verbose {false};
};	

#endif