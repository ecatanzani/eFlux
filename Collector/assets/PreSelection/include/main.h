#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>

#pragma once

struct in_pars
{
	std::string input_path;
	std::string output_path;
	std::string logs_dir;
	bool verbose {false};

	// Tasks
	bool mc {false};
};

#endif