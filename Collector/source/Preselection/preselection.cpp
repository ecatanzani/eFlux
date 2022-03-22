#include "Preselection/preselection.h"

void preselection::Pipeline(
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
    const trigger_info &evt_trigger_info)
    {
        if (evt_trigger_info.HET)
        {
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
                output.maxelayer_cut = true;
            if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                output.maxbarlayer_cut = true;
            if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                output.bgotrack_cut = true;
            
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                    if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                        output.bgofiducial_cut = true;
            
            if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarIndex(13), bgoTotalE))
                output.nbarlayer13_cut = true;

            if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                output.maxrms_cut = true;

            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts)) 
            {
                output.trackselection_cut = true;
                stk_charge_measurement(stkclusters);
                if (psd_fiducial_volume_cut())
                {
                    if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo())) 
                    {
                        output.psdstkmatch_cut = true;
                        psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
                        if (psd_charge_cut(cuts))
                        {
                            output.psdcharge_cut = true;
                            if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                output.stkcharge_cut = true;
                        }            
                    }
                }
                else
                {
                    if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                        output.stkcharge_cut = true;
                } 
            }
            
            reset_cuts_results();

            // maxelayer lastcut
            if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
            {
                if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                {
                    if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                    {
                        if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                        {
                            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                            {
                                stk_charge_measurement(stkclusters);
                                if (psd_fiducial_volume_cut())
                                {
                                    if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                    {
                                        psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
                                        if (psd_charge_cut(cuts))
                                        {
                                            if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                                output.maxelayer_lastcut = true;
                                        }
                                    }
                                }
                                else
                                {
                                    if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                        output.maxelayer_lastcut = true;
                                }
                            }
                        }
                    }
                }
            }

            reset_cuts_results();

            // maxbarlayer lastcut
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
            {
                if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                {
                    if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                    {
                        if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                        {
                            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                            {
                                stk_charge_measurement(stkclusters);
                                if (psd_fiducial_volume_cut())
                                {
                                    if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                    {
                                        psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
                                        if (psd_charge_cut(cuts))
                                        {
                                            if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                                output.maxbarlayer_lastcut = true;
                                        }
                                    }
                                }
                                else
                                {
                                    if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                        output.maxbarlayer_lastcut = true;
                                }
                            }
                        }
                    }
                }
            }

            reset_cuts_results();

            // bgotrack lastcut
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
            {
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                {
                    if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                    {
                        if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                        {
                            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                            {
                                stk_charge_measurement(stkclusters);
                                if (psd_fiducial_volume_cut())
                                {
                                    if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                    {
                                        psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
                                        if (psd_charge_cut(cuts))
                                        {
                                            if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                                output.bgotrack_lastcut = true;
                                        }
                                    }
                                }
                                else
                                {
                                    if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                        output.bgotrack_lastcut = true;
                                }
                            }
                        }
                    }
                }
            }
                                

            reset_cuts_results();

            // BGO fiducial lastcut
            if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
            {
                if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                {
                    if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                    {
                        stk_charge_measurement(stkclusters);
                        if (psd_fiducial_volume_cut())
                        {
                            if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                            {
                                psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
                                if (psd_charge_cut(cuts))
                                {
                                    if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                        output.bgofiducial_lastcut = true;
                                }
                            }
                        }
                        else
                        {
                            if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                output.bgofiducial_lastcut = true;
                        } 
                    }
                }
            }      

            reset_cuts_results();

            // nbarlayer13 lastcut
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
            {
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                {
                    if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                    {
                        if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                        {
                            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                            {
                                stk_charge_measurement(stkclusters);
                                if (psd_fiducial_volume_cut())
                                {
                                    if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                    {
                                        psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
                                        if (psd_charge_cut(cuts))
                                        {
                                            if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                                output.nbarlayer13_lastcut = true;
                                        }
                                    } 
                                }
                                else
                                {
                                    if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                        output.nbarlayer13_lastcut = true;
                                }
                            }
                        }
                    }
                }
            }
                                      
            
            reset_cuts_results();

            // maxrms lastcut
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
            {
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                {
                    if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                    {
                        if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                        {
                            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                            {
                                stk_charge_measurement(stkclusters);
                                if (psd_fiducial_volume_cut())
                                {
                                    if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                    {
                                        psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
                                        if (psd_charge_cut(cuts))
                                        {
                                            if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                                output.maxrms_lastcut = true;
                                        }
                                    }
                                }
                                else
                                {
                                    if (stk_charge_cut(cuts.STK_charge_upper, cuts.STK_charge_medium, cuts.STK_charge_lower))
                                        output.maxrms_lastcut = true;
                                }
                            }
                        }
                    }
                }
            }
            
            reset_cuts_results();

            // trackselection lastcut
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                    if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                        if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                            if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                                output.trackselection_lastcut = true;

            reset_cuts_results();

            // psdstkmatch lastcut
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                    if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                        if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                            if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                                if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                                    output.psdstkmatch_lastcut = true;
        }
    }

const presel_output preselection::GetPreselectionOutput()
{
    return output;
}

void preselection::reset_preselection_output()
{
    output.maxelayer_cut            = false;
    output.maxbarlayer_cut          = false;
    output.bgotrack_cut             = false;
    output.bgofiducial_cut          = false;
    output.nbarlayer13_cut          = false;
    output.maxrms_cut               = false;
    output.trackselection_cut       = false;
    output.psdstkmatch_cut          = false;
    output.psdcharge_cut            = false;
    output.stkcharge_cut            = false;
    output.maxelayer_lastcut        = false;
    output.maxbarlayer_lastcut      = false;
    output.bgotrack_lastcut         = false;
    output.bgofiducial_lastcut      = false;
    output.nbarlayer13_lastcut      = false;
    output.maxrms_lastcut           = false;
    output.trackselection_lastcut   = false;
    output.psdstkmatch_lastcut      = false;
}

void preselection::reset_cuts_results()
{
    reset_stk_best_track();
	reset_psd_clusters();
	reset_psd_charges();
	reset_stk_charges();
	reset_classifiers();
	reset_filter_output();
	reset_time();
	reset_trigger();
}

void preselection::Reset()
{
    reset_cuts_results();
    reset_preselection_output();
}