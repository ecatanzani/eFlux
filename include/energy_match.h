#ifndef ENERGY_MATCH_H
#define ENERGY_MATCH_H

#include <iostream>
#include <vector>

extern void allocateParticleEnergy(
    std::vector<unsigned int> &counts,
    const std::vector<float> &logEBins,
    const double totalE);

extern int linearSearch(
    const std::vector<float> &logEBins, 
    const double energy);

extern int binarySearch(
    const std::vector<float> &logEBins, 
    int l, 
    int r, 
    const double energy);

#endif