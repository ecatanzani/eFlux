#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <memory>

#include "anyoption.h"

#pragma once

struct in_args
{
    std::string collector_wd;
	std::string input_list;
	std::string output_path;
	std::string fit_tree_path;
	bool verbose {false};
	bool mc_flag {false};
	unsigned int threads {1};
};

#endif