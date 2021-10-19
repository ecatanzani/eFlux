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
	std::string config_dir;
	std::string output_path;
	std::string learning_method;
	bool verbose = false;
};

#endif