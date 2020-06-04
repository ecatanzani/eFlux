#ifndef MYHEADER_H
#define MYHEADER_H

#include <iostream>
#include <string>
#include <vector>

#include "anyoption.h"

#pragma once

extern void eCore(
    const std::string inputPath,
    const std::string outputPath,
    const bool verbose,
    const bool pedantic,
    const unsigned int lvTime,
    const std::string accInputPath,
    AnyOption &opt,
    const std::string wd);

extern const std::string uniqueOutFile(
    const std::string outputPath,
    AnyOption &opt);

extern bool chechFlags(
    AnyOption &opt,
    const std::string inputPath,
    const std::string outputPath,
    const unsigned int lvTime);

extern void generateFinalGraph(
    const bool verbose,
    const bool pedantic,
    const std::string outputPath,
    const std::string complete_histo_path,
    const std::string wd);

#endif