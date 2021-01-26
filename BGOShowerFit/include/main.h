#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <memory>

#include "anyoption.h"

#pragma once

void reader(
    const std::string wd,
    const std::string inputList,
    const std::string outputPath,
    const bool _VERBOSE,
    const bool mc);

#endif