#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include "Dmp/DmpFilterContainer.h"

struct eff_output 
{
    bool trigger_efficiency_preselection                                    {false};
    bool trigger_efficiency_preselection_is_het                             {false};
    bool trigger_efficiency_preselection_is_let                             {false};
    bool trigger_efficiency_preselection_is_unb                             {false};

    bool maxrms_efficiency_preselection                                     {false};
    bool maxrms_efficiency_preselection_accepted                            {false};

    bool nbarlayer13_efficiency_preselection                                {false};
    bool nbarlayer13_efficiency_preselection_accepted                       {false};

    bool maxrms_and_nbarlayer13_efficiency_preselection                     {false};
    bool maxrms_and_nbarlayer13_efficiency_preselection_accepted            {false};

    bool sumrms_low_energy_cut_efficiency_preselection                      {false};
    bool sumrms_low_energy_cut_efficiency_preselection_accepted             {false};

    bool track_efficiency_preselection                                      {false};
    bool track_efficiency_preselection_accepted                             {false};

    bool stk_1rm_cut_efficiency_preselection                                {false};
    bool stk_1rm_cut_efficiency_preselection_accepted                       {false};

    bool psdstkmatch_efficiency_preselection                                {false};
    bool psdstkmatch_efficiency_preselection_accepted                       {false};

    bool psdcharge_efficiency_preselection                                  {false};
    bool psdcharge_efficiency_preselection_accepted                         {false};

    bool stkcharge_efficiency_preselection                                  {false};
    bool stkcharge_efficiency_preselection_accepted                         {false};
};

class efficiency: public DmpFilterContainer
{
    public:
        efficiency(){};
        ~efficiency(){};

        void Pipeline(
            const std::shared_ptr<DmpEvtBgoRec> &bgorec,
            const std::shared_ptr<DmpEvtBgoHits> &bgohits,
            const cuts_conf &cuts,
            const double bgoTotalE,
            const double bgoTotalE_corr,
            DmpBgoContainer &bgoVault,
            DmpStkContainer &stkVault,
            DmpPsdContainer &psdVault,
            const std::shared_ptr<TClonesArray> &stkclusters,
            const std::shared_ptr<TClonesArray> &stktracks,
            const active_cuts &acuts,
            const trigger_info &evt_trigger_info);

        const eff_output GetEfficiencyOutput();
        void Reset();
    
    private:
        void reset_cuts_results();
        void reset_efficiency_output();
        eff_output output;
        
};

#endif