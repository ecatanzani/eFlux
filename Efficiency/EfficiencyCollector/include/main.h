#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>

#pragma once

struct in_pars
{
    std::string wd;
	std::string input_path;
	std::string output_path;

	bool verbose = false;
    unsigned int threads = 1;

    // Tasks
	bool mc_flag = false;
	bool rawdata_flag = false;

    bool CheckArgs()
    {
        bool status = true;
        if (wd.empty() || input_path.empty() || output_path.empty())
        {
            status = false;
            std::cerr << "\n(ERROR) -> Please check input parameters, they might be empty !\n";
        }
        if (mc_flag==false && rawdata_flag==false)
        {
            status = false;
            std::cerr << "\n(ERROR) -> Please choose a task... \n";
        }
        if (mc_flag==true && rawdata_flag==true)
        {
            status = false;
            std::cerr << "\n(ERROR) -> Only one task can be choosen... \n";
        }
        if (!status)
            std::cout << std::endl;
        return status;
    }

    void ExpandArgs()
    {
        std::cout << "\n\n***** Input Args *****\n\n";
        std::cout << "Input file list [" << input_path << "]";
        std::cout << "\nCollector config directory [" << wd << "]";
        std::cout << "\nOutput path [" << output_path << "]";
        verbose ? std::cout << "\nVerbosity [True]" : std::cout << "\nVerbosity [False]";
        std::cout << "\n\n**********************\n";
    }
};

#endif