#ifndef TIMING_H
#define TIMING_H

#include <map>
#include <utility>

struct timing
{
    unsigned int prev_sec = 0; 
    unsigned int curr_sec = 0;
    std::map<unsigned int, bool> repeated_sec;

    void SetPrevSec(const unsigned int sec)
    {
        prev_sec = sec;
    }
    void SetCurrSec(const unsigned int sec)
    {
        curr_sec = sec;
    }
    void CheckRepeated()
    {
        if (curr_sec<prev_sec)
        {
            repeated_sec.insert(std::make_pair(prev_sec,false));
            repeated_sec.insert(std::make_pair(curr_sec,false));
        }        
    }
};

#endif