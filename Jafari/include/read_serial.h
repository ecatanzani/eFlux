#ifndef READ_SERIAL_H
#define READ_SERIAL_H

#include <string>
#include <iostream>

#include "TFile.h"

extern void read_serial(
    const std::string inputPath,
    TFile& outFile,
    const std::string wd,
    const bool verbose);

#endif