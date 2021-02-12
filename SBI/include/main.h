#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>

#include "anyoption.h"

#pragma once

struct in_pars
{
	std::string input_path;
	std::string output_path;

	bool verbose = false;
	bool pedantic = false;

    bool CheckArgs()
    {
        bool status = true;
        if (input_path.empty() || output_path.empty())
        {
            status = false;
            std::cerr << "\n(ERROR) -> Please check input parameters, they might be empty !\n";
        }
        return status;
    }

    void ExpandArgs()
    {
        std::cout << "\n\n***** Input Args *****\n\n";
        std::cout << "Input file list [" << input_path << "]";
        std::cout << "\nOutput path [" << output_path << "]";
        verbose ? std::cout << "\nVerbosity [True]" : std::cout << "\nVerbosity [False]";
        pedantic ? std::cout << "\nPedantic [True]" : std::cout << "\nPedantic [False]";
        std::cout << "\n\n**********************\n";
    }
};

#endif