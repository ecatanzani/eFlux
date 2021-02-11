#ifndef CONTAINER_H
#define CONTAINER_H

#include <map>
#include <vector>
#include <utility>
#include <iostream>

#include "DmpEvtNudRaw.h"
#include "DmpEvtHeader.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtPsdHits.h"
#include "DmpEvtStkReco.hh"
#include "DmpEvtAttitude.h"

#include "TClonesArray.h"

#define _trig 5
#define _msdeadtime 3.0725

struct stats
{
    unsigned int processed_events = 0;
    unsigned int total_events;
    unsigned int seconds = 0;
    unsigned int good_sbi_events = 0;
    unsigned int bad_sbi_events = 0;
};

struct timing
{
    unsigned int prev_sec = 0;
    unsigned int curr_sec = 0;
    bool newsec = false;

    std::map<unsigned int, bool> repeated_sec;

    // Set the previous second
    void SetPrevSec(const unsigned int sec)
    {
        prev_sec = sec;
    }

    /*
        - Set current second 
        - Checks if the second cnahges respect to the previous one
        - Checks if the seconds has been repeated
    */
    void SetCurrSec(const unsigned int sec)
    {   
        curr_sec = sec;
        if (curr_sec != prev_sec)
            newsec = true;
        else
            newsec = false;
        CheckRepeated();
        
    }
    void CheckRepeated()
    {
        if (curr_sec < prev_sec)
        {
            repeated_sec.insert(std::make_pair(prev_sec, false));
            repeated_sec.insert(std::make_pair(curr_sec, false));
        }
        prev_sec = curr_sec;
    }
    bool CheckNewSec()
    {
        return newsec;
    }
    void ShowSecs()
    {
        std::cout << "\n\nPrevous sec: " << prev_sec;
        std::cout << "\nCurrent sec: " << curr_sec;
    }
};

struct container
{
    timing evt_time;
    stats evt_stat;

    bool goodsbi = false;
    unsigned int run = 0;
    unsigned int second = 0;
    double lvtime = 0;
    unsigned int nev = 0;
    unsigned int firstev = 0;
    unsigned int lastev = 0;
    std::vector<unsigned int> ntrig;
    std::vector<bool> trigenable;
    unsigned int ntrigperiod = 0;
    bool trigperiodenable = false;
    bool stkstatus = false;
    bool bgostatus = false;
    bool psdstatus = false;
    bool nudstatus = false;
    unsigned int ntracks = 0;
    unsigned int nbgohits = 0;
    unsigned int npsdhits = 0;
    unsigned int nnudhits = 0;

    double msdtime = 0; // deadtime to perform operations
    double msdtime_beginning = 0; //dead time to add at the beginning od the second, to recover events triggered at the end of the previous second

    bool insaa = false;

    container() :   ntrig (_trig, 0),
                    trigenable (_trig, false)
    {
    }

    bool SetSBIStatus(
        std::shared_ptr<DmpEvtHeader> header,
        std::shared_ptr<DmpEvtAttitude> attitude,
        const unsigned int idx)
    {
        bool status = true;

        // Header check
        if (header == nullptr)
        {
            std::cerr << "\nWARNING: Header obj not found at entry [" << idx << "] -- skipping\n";
            goodsbi = false;
            status = false;
        }
        else
        {
            if (idx!=1)
                goodsbi = true;
        }

        // Attitude check
        if (attitude == nullptr)
        {
            std::cerr << "\nWARNING: Attitude obj not found at entry [" << idx << "] -- skipping\n";
            goodsbi = false;
        }

        return status;
    }

    void SetSecond()
    {
        second = evt_time.curr_sec;
    }
    void SetEventNumber(const unsigned int ev)
    {
        nev = ev;
    }
    void SetLastEventNumber(const unsigned int ev)
    {
        lastev = ev;
    }
    void UpdateEventNumber()
    {
        ++nev;
    }
    void UpdateRunNumber()
    {
        run = second;
    }
    void UpdateDeadTime(const unsigned int msecond)
    {
        msdtime += _msdeadtime;
        if ((double)msecond>(1000-_msdeadtime))
        {
            msdtime -= (_msdeadtime -(1000-(double)msecond));
            msdtime_beginning += (_msdeadtime -(1000-(double)msecond) );
        }
    }
    void UpdateTrigger(std::shared_ptr<DmpEvtHeader> header)
    {
        for (int itrig=0; itrig<_trig; itrig++)
        {
            if (header->GeneratedTrigger(itrig))
                ++ntrig[itrig];
            if (!header->EnabledTrigger(itrig))
                trigenable[itrig]=false;
        }
        
        if (header->GeneratedPeriodTrigger())
            ++ntrigperiod;
        if (!header->EnabledPeriodTrigger())
            trigperiodenable=false;
    }
    void UpdateSubDetectorsOccupancy(
        std::shared_ptr<TClonesArray> stktracks,
        std::shared_ptr<DmpEvtBgoRec> bgorec,
        std::shared_ptr<DmpEvtPsdHits> psdhits,
        std::shared_ptr<DmpEvtNudRaw> nudraw)
    {
        auto nudtotalhits = [](short ev_ch0, short ev_ch1, short ev_ch2, short ev_ch3) -> short { return ev_ch0 + ev_ch1 + ev_ch2 + ev_ch3; };

        ntracks  += stktracks->GetLast() + 1;
        nbgohits += bgorec->GetTotalHits();
        npsdhits += psdhits->GetHittedBarNumber();
        nnudhits += (unsigned int)nudtotalhits(nudraw->fChannelID[0], nudraw->fChannelID[1], nudraw->fChannelID[2], nudraw->fChannelID[3]);
    }
    void UpdateSubDetectorStatus(std::shared_ptr<DmpEvtHeader> header)
    {
        if (header->GetSubDetectorStatus(DmpEDetectorID::kStk)) 
            stkstatus=false;
        if (header->GetSubDetectorStatus(DmpEDetectorID::kBgo)) 
            bgostatus=false;
        if (header->GetSubDetectorStatus(DmpEDetectorID::kPsd))
            psdstatus=false;
        if (header->GetSubDetectorStatus(DmpEDetectorID::kNud))
            nudstatus=false;
    }
    void CheckSAA(std::shared_ptr<DmpEvtAttitude> attitude)
    {
        insaa |= (attitude!=nullptr) ? (bool)(attitude->insaa) : 0;
    }
    bool CheckRepeatedSeconds()
    {
        bool fill = true;
        auto it = evt_time.repeated_sec.find(second);
        if (it != evt_time.repeated_sec.end())
        {
            if (it->second==false)
            {
	            it->second=true;
	            goodsbi &= false;
            }
            else
	            fill=false;
        }
        return fill;
    }
    void ComputeLiveTime()
    {
        lvtime = (1000 - msdtime)/(double)1000;
    }
    void CleanUp(const unsigned int first_ev)
    {
        msdtime = msdtime_beginning;
        msdtime_beginning = 0;
        Reset(first_ev);
    }
    void Reset(const unsigned int first_ev)
    {
        goodsbi = false;
        run = 0;
        second = evt_time.curr_sec;
        lvtime = 0;
        nev = 0;
        firstev = first_ev;
        lastev = 0;
        ntrig = std::vector<unsigned int> (_trig, 0);
        trigenable = std::vector<bool> (_trig, 0);
        ntrigperiod = 0;
        trigperiodenable = false;
        stkstatus = false;
        bgostatus = false;
        psdstatus = false;
        nudstatus = false;
        ntracks = false;
        nbgohits = 0;
        npsdhits = 0;
        nnudhits = 0;
    }
    void SetSBIStatus(const bool status)
    {
        goodsbi = status;
    }
    void UpdateStats(const bool good_evt=true)
    {
        if (good_evt)
            ++evt_stat.processed_events;
        ++evt_stat.total_events;
        if (goodsbi)
            ++evt_stat.good_sbi_events;
        else
            ++evt_stat.bad_sbi_events;
    }
    void UpdateNewSecStats()
    {
        ++evt_stat.seconds;
    }
    void PrintStats()
    {
        std::cout << "\n\n *** Stats ***\n\n";
        std::cout << "Total number of events: " << evt_stat.total_events << std::endl;
        std::cout << "Total number of processed events: " << evt_stat.processed_events << std::endl;
        std::cout << "Total number of GOOD sbi events: " << evt_stat.good_sbi_events << std::endl;
        std::cout << "Total number of BAD sbi events: " << evt_stat.bad_sbi_events << std::endl;
        std::cout << "Total number of seconds: " << evt_stat.seconds << std::endl;
        std::cout << "\n**********\n";
    }
};

#endif