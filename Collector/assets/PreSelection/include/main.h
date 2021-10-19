#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>

#pragma once

struct in_pars
{
	std::string input_path;
	std::string output_path;
	bool verbose {false};
    unsigned int threads {1};

	// Tasks
	bool mc {false};
};

#endif