#include "Efficiency/efficiency.h"

void efficiency::Pipeline(
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
        // Trigger efficiency
        if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
            if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                    if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                        if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                                if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                    {
                                        psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
                                        if (psd_charge_cut(cuts))
                                        {
                                            output.trigger_efficiency_preselection = true;
                                            if (evt_trigger_info.HET) output.trigger_efficiency_preselection_is_het = true;
                                            if (evt_trigger_info.LET) output.trigger_efficiency_preselection_is_let = true;
                                            if (evt_trigger_info.unbiased) output.trigger_efficiency_preselection_is_unb = true;
                                        }
                                    }
        
        reset_cuts_results();
        
        // MaxRms efficiency
        if (evt_trigger_info.HET)
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                    if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                        if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                                if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                {
                                    psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
                                    if (psd_charge_cut(cuts))
                                    {
                                        output.maxrms_efficiency_preselection = true;
                                        if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                                            output.maxrms_efficiency_preselection_accepted = true;
                                    }
                                }
        
        reset_cuts_results();

        // nbarlayer13 efficiency
        if (evt_trigger_info.HET)
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                    if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                        if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                                if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                {
                                    psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
                                    if (psd_charge_cut(cuts))
                                    {
                                        output.nbarlayer13_efficiency_preselection = true;
                                        if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                                            output.nbarlayer13_efficiency_preselection_accepted = true;
                                    }
                                }

        reset_cuts_results();

        // MaxRms && nbarlayer13 efficiency
        if (evt_trigger_info.HET)
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                    if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                        if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                            if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                            {
                                psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
                                if (psd_charge_cut(cuts))
                                {
                                    output.maxrms_and_nbarlayer13_efficiency_preselection = true;
                                    if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts) && nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                                        output.maxrms_and_nbarlayer13_efficiency_preselection_accepted = true;
                                }
                            }

        reset_cuts_results();

        // track selection efficiency
        if (evt_trigger_info.HET)
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                    if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                        if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                            if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                            {
                                output.track_efficiency_preselection = true;
                                if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                                    output.track_efficiency_preselection_accepted = true;
                            }

        reset_cuts_results();

        // psd_stk_match efficiency
        if (evt_trigger_info.HET)
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                    if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                        if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                            if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                                if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                                {
                                    output.psdstkmatch_efficiency_preselection = true;
                                    if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                        output.psdstkmatch_efficiency_preselection_accepted = true;
                                }

        reset_cuts_results();

        // psd_charge efficiency
        if (evt_trigger_info.HET)
            if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
                if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                    if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                        if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                            if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                                if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                                    if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                    {
                                        psd_charge_measurement(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID());
						                stk_charge_measurement(stkclusters);
                                        if (stk_charge_cut(11))
                                        {
                                            output.psdcharge_efficiency_preselection = true;
                                            if (psd_charge_cut(cuts))
                                                output.psdcharge_efficiency_preselection_accepted = true;
                                        }
                                    }
    }

const eff_output efficiency::GetEfficiencyOutput()
{
    return output;
}

void efficiency::reset_efficiency_output()
{
    output.trigger_efficiency_preselection                                  = false;
    output.trigger_efficiency_preselection_is_het                           = false;
    output.trigger_efficiency_preselection_is_let                           = false;
    output.trigger_efficiency_preselection_is_unb                           = false;

    output.maxrms_efficiency_preselection                                   = false;
    output.maxrms_efficiency_preselection_accepted                          = false;

    output.nbarlayer13_efficiency_preselection                              = false;
    output.nbarlayer13_efficiency_preselection_accepted                     = false;

    output.maxrms_and_nbarlayer13_efficiency_preselection                   = false;
    output.maxrms_and_nbarlayer13_efficiency_preselection_accepted          = false;

    output.track_efficiency_preselection                                    = false;
    output.track_efficiency_preselection_accepted                           = false;

    output.psdstkmatch_efficiency_preselection                              = false;
    output.psdstkmatch_efficiency_preselection_accepted                     = false;

    output.psdcharge_efficiency_preselection                                = false;
    output.psdcharge_efficiency_preselection_accepted                       = false;
}


void efficiency::reset_cuts_results()
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

void efficiency::Reset()
{
    reset_cuts_results();
    reset_efficiency_output();
}