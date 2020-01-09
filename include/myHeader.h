#ifndef MYHEADER_H
#define MYHEADER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TTree.h"

#pragma once

extern void eCore(
                    const std::string inputPath,
                    const std::string outputPath,
                    const bool verbose
                );

extern void readInputTree(const std::string inputPath,std::vector<double> &dataValues);
extern void branchTree(TTree &myDataTree,std::vector<double> &dataValues);

#endif