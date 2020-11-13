#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <memory>

#include "config.h"
#include "anyoption.h"

#include "TChain.h"

#pragma once

void reader(
    const std::string wd,
    const std::string inputList,
    const std::string outputPath,
    const bool _VERBOSE,
    const bool mc);

void mc_reader(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE);

#endif