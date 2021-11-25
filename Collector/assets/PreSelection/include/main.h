#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>

#include "anyoption.h"

#pragma once

struct in_pars
{
	std::string input_path;
	std::string output_wd;
	std::string config_wd;
	bool verbose {false};

	// Tasks
	bool mc_flag {false};
	bool rawdata_flag {false};
	bool rank_tracks {false};
};

#endif