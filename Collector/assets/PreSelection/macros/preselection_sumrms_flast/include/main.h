#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>

#pragma once

struct in_pars
{
	std::string input_path;
    std::string output_directory;
    std::string energy_config_file;
	std::string output_path;
    double input_spectral_index {0};

	bool verbose {false};

	// Tasks
	bool mc_flag {false};

    unsigned int threads {1};

    bool CheckArgs()
    {
        bool status = true;
        if (input_path.empty() || output_path.empty())
        {
            status = false;
            std::cerr << "\n(ERROR) -> Please check input parameters, they might be empty !\n";
        }
        if (!status)
            std::cout << std::endl;
        return status;
    }

    void ExpandArgs()
    {
        std::cout << "\n\n***** Input Args *****\n\n";
        std::cout << "Input file list [" << input_path << "]";
        std::cout << "\nOutput path [" << output_path << "]";
        verbose ? std::cout << "\nVerbosity [True]" : std::cout << "\nVerbosity [False]";
        mc_flag ? std::cout << "\nMC flag [True]" : std::cout << "\nMC flag [False]";
        std::cout << "\n\n**********************\n";
    }

    void GetOutputFile()
    {
        output_path = output_directory + std::string("/preselection.root");
    }
};

#endif