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
		const active_cuts &acuts)
    {
        // Trigger efficiency
        if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
            if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                    if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                        if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                                if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
							        if (psd_charge_cut(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID(), cuts))
                                    {
                                        output.trigger_efficiency_preselection = true;
                                        if (evt_trigger.HET) output.trigger_efficiency_preselection_is_het = true;
                                        if (evt_trigger.LET) output.trigger_efficiency_preselection_is_let = true;
                                        if (evt_trigger.unbiased) output.trigger_efficiency_preselection_is_unb = true;
                                    }

        // MaxRms efficiency
        if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
            if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                    if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                        if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                            if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                if (psd_charge_cut(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID(), cuts))
                                {
                                    output.maxrms_efficiency_preselection = true;
                                    if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                                        output.maxrms_efficiency_preselection_accepted = true;
                                    else
                                        output.maxrms_efficiency_preselection_notaccepted = true;
                                }

        // nbarlayer13 efficiency
        if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
            if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                    if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                        if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                            if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                if (psd_charge_cut(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID(), cuts))
                                {
                                    output.nbarlayer13_efficiency_preselection = true;
                                    if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                                        output.nbarlayer13_efficiency_preselection_accepted = true;
                                    else
                                        output.nbarlayer13_efficiency_preselection_notaccepted = true;
                                }

        // track selection efficiency
        if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
            if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                    if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                        if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                        {
                            output.track_efficiency_preselection = true;
                            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                                output.track_efficiency_preselection_accepted = true;
                            else
                                output.track_efficiency_preselection_notaccepted = true;
                        }
                                        
        // psd_stk_match efficiency
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
                                else
                                    output.psdstkmatch_efficiency_preselection_notaccepted = true;
                            }

        // psd_charge efficiency
        if (maxElayer_cut(bgoVault.GetLayerEnergies(), cuts, bgoTotalE))
            if (maxBarLayer_cut(bgoVault.GetLayerBarNumber(), bgoVault.GetiMaxLayer(), bgoVault.GetIdxBarMaxLayer()))
                if (BGOTrackContainment_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts))
                    if (nBarLayer13_cut(bgohits, bgoVault.GetSingleLayerBarNumber(13), bgoTotalE))
                        if (maxRms_cut(bgoVault.GetELayer(), bgoVault.GetRmsLayer(), bgoTotalE, cuts))
                            if (track_selection_cut(bgorec, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), bgohits, stkclusters, stktracks, cuts))
                                if (psd_stk_match_cut(bgoVault.GetBGOslope(), bgoVault.GetBGOintercept(), cuts, psdVault.getPsdClusterIdxBegin(), psdVault.getPsdClusterZ(), psdVault.getPsdClusterMaxECoo()))
                                {
                                    output.psdcharge_efficiency_preselection = true;
                                    if (psd_charge_cut(psdVault.getPsdClusterMaxE(), psdVault.getPsdClusterIdxMaxE(), psdVault.getHitZ(), psdVault.getGlobalBarID(), cuts))
                                        output.psdcharge_efficiency_preselection_accepted = true;
                                    else
                                        output.psdcharge_efficiency_preselection_notaccepted = true;
                                }
    }

const eff_output efficiency::GetEfficiencyOutput()
{
    return output;
}

void efficiency::reset_efficiency_output()
{
    output.trigger_efficiency_preselection                       = false;
    output.trigger_efficiency_preselection_is_het                = false;
    output.trigger_efficiency_preselection_is_let                = false;
    output.trigger_efficiency_preselection_is_unb                = false;

    output.maxrms_efficiency_preselection                        = false;
    output.maxrms_efficiency_preselection_accepted               = false;
    output.maxrms_efficiency_preselection_notaccepted            = false;

    output.nbarlayer13_efficiency_preselection                   = false;
    output.nbarlayer13_efficiency_preselection_accepted          = false;
    output.nbarlayer13_efficiency_preselection_notaccepted       = false;

    output.track_efficiency_preselection                         = false;
    output.track_efficiency_preselection_accepted                = false;
    output.track_efficiency_preselection_notaccepted             = false;

    output.psdstkmatch_efficiency_preselection                   = false;
    output.psdstkmatch_efficiency_preselection_accepted          = false;
    output.psdstkmatch_efficiency_preselection_notaccepted       = false;

    output.psdcharge_efficiency_preselection                     = false;
    output.psdcharge_efficiency_preselection_accepted            = false;
    output.psdcharge_efficiency_preselection_notaccepted         = false;
}

void efficiency::Reset()
{
    reset_stk_best_track();
	reset_psd_clusters();
	reset_psd_charges();
	reset_stk_charges();
	reset_classifiers();
	reset_filter_output();
	reset_time();
	reset_trigger();
    reset_efficiency_output();
}