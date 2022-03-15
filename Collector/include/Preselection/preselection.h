#ifndef PRESELECTION_H
#define PRESELECTION_H

#include "Dmp/DmpFilterContainer.h"

struct presel_output
{
    // Cut variables
    bool maxelayer_cut              {false};
    bool maxbarlayer_cut            {false};
    bool bgotrack_cut               {false};
    bool bgofiducial_cut            {false};
    bool nbarlayer13_cut            {false};
    bool maxrms_cut                 {false};
    bool trackselection_cut         {false};
    bool psdstkmatch_cut            {false};
    bool psdcharge_cut              {false};
    bool stkcharge_cut              {false};
    
    // Last Cut variables
    bool maxelayer_lastcut          {false};        // All cuts exect maxelayer
    bool maxbarlayer_lastcut        {false};        // All cuts except maxbarlayer
    bool bgotrack_lastcut           {false};        // All cuts except bgotrack
    bool bgofiducial_lastcut        {false};        // All cuts except BGO fiducial ones
    bool nbarlayer13_lastcut        {false};        // All cuts except nbarlayer13
    bool maxrms_lastcut             {false};        // All cuts except maxrms
    bool trackselection_lastcut     {false};        // All cuts except trackselection (and the related ones)
    bool psdstkmatch_lastcut        {false};        // All cuts except psd-stk matching (and the related ones)
};

class preselection: public DmpFilterContainer
{
    public:
        preselection(){};
        ~preselection(){};

        void Pipeline(
            const std::shared_ptr<DmpEvtBgoRec> &bgorec,
            const std::shared_ptr<DmpEvtBgoHits> &bgohits,
            const cuts_conf &cuts,
            const double bgoTotalE,
            const double bgoTotalE_corr,
            DmpBgoContainer &bgoVault,
            DmpPsdContainer &psdVault,
            const std::shared_ptr<TClonesArray> &stkclusters,
            const std::shared_ptr<TClonesArray> &stktracks,
            const active_cuts &acuts,
            const trigger_info &evt_trigger_info);
        
        const presel_output GetPreselectionOutput();
        void Reset();

    private:
        void reset_cuts_results();
        void reset_preselection_output();
        presel_output output;
};

#endif