#include "energy_match.h"
#include "read_sets_config_file.h"

void allocateParticleEnergy(
    std::vector<unsigned int> &counts,
    const std::vector<float> &logEBins,
    const double totalE)
{
    const double energy = totalE / 1000;                               // Scale energy in GeV - binning is in GeV !!
    auto idx = binarySearch(logEBins, 0, logEBins.size() - 1, energy); //More efficient vector search respect to the linear one
    //auto idx = linearSearch(logEBins,energy);
    if (idx != -1)
    {
        /*
        if ( !(energy > logEBins[idx] && energy < logEBins[idx+1]) ) 
            std:: cout << "\n Idx: " << idx << "\tLow(" << logEBins[idx] << "): " << idx << "\tHigh(" << logEBins[idx+1] << "): " << idx+1 << "\tEgy: " << energy;
        */
        ++counts[idx];
    }
    else
        std::cerr << "\nWARNING: Could not find energy bin for current event: " << energy << " GeV";

    /*
        USE EXPECTIONS FOR THAT
    */
}

int linearSearch(
    const std::vector<float> &logEBins, 
    const double energy)
{
    /*
        ****** ALL VECTOR SCAN ******
    */

    unsigned int lIdx = 0;
    bool found_interval = false;
    for (unsigned int vIdx = 0; vIdx < logEBins.size() - 1; ++vIdx)
        if (energy > logEBins[vIdx] && energy < logEBins[vIdx + 1])
        {
            lIdx = vIdx;
            found_interval = true;
            break;
        }
    //std::cout << "\nEnergy: " << energy << "\tLow: " << logEBins[idx] << "\tHigh: " << logEBins[idx+1];
    if (found_interval)
        return lIdx;
    else
        return -1;
}

int binarySearch(const std::vector<float> &logEBins, int l, int r, const double energy)
{
    /*
        ****** BINARY SEARCH ******
    */

    if (r >= l)
    {
        int mid = l + (r - l) / 2;
        if (energy > logEBins[mid] && energy < logEBins[mid + 1])
            return mid;
        if (logEBins[mid] > energy)
            return binarySearch(logEBins, l, mid - 1, energy);
        return binarySearch(logEBins, mid + 1, r, energy);
    }
    return -1;
}

int getInputPowerLawIndex(const double lEnergy, const double hEnergy, const data_set_conf input_sets)
{
    int plawIndex = 0;
    for(auto it=input_sets.sets_eMin.begin(); it!=input_sets.sets_eMin.end(); ++it)
    {
        auto index = std::distance(input_sets.sets_eMin.begin(), it);
        auto tmp_eMin = *it;
        auto tmp_eMax = input_sets.sets_eMax[index];
        if(lEnergy>=tmp_eMin && hEnergy<=tmp_eMax)
        {
            plawIndex = input_sets.sets_powerlaw_idx[index];
            break;
        }
    }
    return plawIndex;
}