#include "acceptance.h"
#include "acceptance_cuts.h"
#include "data_cuts.h"

bool filter_this_mc_event(
    event_filter &filter,
    const std::shared_ptr<DmpEvtSimuPrimaries> &simu_primaries,
    const std::shared_ptr<DmpEvtBgoRec> &bgorec,
    const std::shared_ptr<DmpEvtBgoHits> &bgohits,
    const cuts_conf &acceptance_cuts,
    const double bgoTotalE,
	DmpBgoContainer &bgoVault,
	DmpPsdContainer &psdVault,
	psd_charge &extracted_psd_charge,
	stk_charge &extracted_stk_charge,
    const std::shared_ptr<TClonesArray> &stkclusters,
    const std::shared_ptr<TClonesArray> &stktracks,
	const data_active_cuts &active_cuts)
{	
	// Create best_track struct
	best_track event_best_track;
		
	// Create PSD-STK match struct
	psd_cluster_match clu_matching;
		
	// **** Geometric cut ****
	filter.geometric = geometric_cut(simu_primaries);

	// **** BGO Fiducial Volume ****
	// maxElayer_cut
	filter.BGO_fiducial_maxElayer_cut = maxElayer_cut(
		bgorec, 
		acceptance_cuts, 
		bgoTotalE);
	// maxBarLayer_cut
	filter.BGO_fiducial_maxBarLayer_cut = maxBarLayer_cut(
		bgoVault.GetLayerBarNumber(), 
		bgoVault.GetiMaxLayer(), 
		bgoVault.GetIdxBarMaxLayer());
	// BGOTrackContainment_cut
	filter.BGO_fiducial_BGOTrackContainment_cut = BGOTrackContainment_cut(
		bgorec, 
		acceptance_cuts);
	
	filter.BGO_fiducial = filter.BGO_fiducial_maxElayer_cut && filter.BGO_fiducial_maxBarLayer_cut && filter.BGO_fiducial_BGOTrackContainment_cut;
	filter.all_cut = filter.BGO_fiducial;

	if (filter.all_cut)
	{
		// **** nBarLayer13 cut ****
		if (active_cuts.nBarLayer13)
		{
			filter.nBarLayer13_cut = nBarLayer13_cut(
				bgohits,
				bgoVault.GetSingleLayerBarNumber(13),
				bgoTotalE);
			filter.all_cut *= filter.nBarLayer13_cut;
		}

		if (filter.all_cut)
		{
			// **** maxRms cut ****
			if (active_cuts.maxRms)
			{
				filter.maxRms_cut = maxRms_cut(
					bgoVault.GetLayerBarNumber(),
					bgoVault.GetRmsLayer(),
					bgoTotalE,
					acceptance_cuts);
				filter.all_cut *= filter.maxRms_cut;
			}

			if (filter.all_cut)
			{
				// **** track selection cut ****
				if (active_cuts.track_selection)
				{
					filter.track_selection_cut = track_selection_cut(
						bgorec,
						bgohits,
						stkclusters,
						stktracks,
						acceptance_cuts,
						event_best_track);
					filter.all_cut *= filter.track_selection_cut;
				}

				if (filter.all_cut)
				{
					// **** psd stk match cut ****
					if (active_cuts.psd_stk_match)
					{
						filter.psd_stk_match_cut = psd_stk_match_cut(
							bgorec,
							acceptance_cuts,
							event_best_track,
							clu_matching,
							psdVault.getPsdClusterIdxBegin(),
							psdVault.getPsdClusterZ(),
							psdVault.getPsdClusterMaxECoo());
						filter.all_cut *= filter.psd_stk_match_cut;
					}

					if (filter.all_cut)
					{	
						// **** psd charge cut ****
						if (active_cuts.psd_charge)
						{
							filter.psd_charge_measurement = true;
							filter.psd_charge_cut = psd_charge_cut(
								event_best_track,
								clu_matching,
								psdVault.getPsdClusterMaxE(),
								psdVault.getPsdClusterIdxMaxE(),
								psdVault.getHitZ(),
								psdVault.getGlobalBarID(),
								acceptance_cuts,
								extracted_psd_charge);
							filter.all_cut *= filter.psd_charge_cut;
						}
						
					
						if (filter.all_cut)
						{
							// **** stk charge cut ****
							if (active_cuts.stk_charge)
							{
								filter.stk_charge_measurement = true;
								filter.stk_charge_cut = stk_charge_cut(
									event_best_track,
									stkclusters,
									acceptance_cuts,
									extracted_stk_charge);
								filter.all_cut *= filter.stk_charge_cut;
							}

							if (filter.all_cut)
							{
								filter.all_cut_no_xtrl = true;

								// **** xtrl cut ****
								if (active_cuts.xtrl)
								{	
									filter.xtrl_cut = xtrl_cut(
										bgoVault.GetSumRMS(),
										bgoVault.GetLastFFracLayer(),
										acceptance_cuts);
									filter.all_cut *= filter.xtrl_cut;
								}
							}
						}
					}
				}		
			}
		}
	}

	if (filter.all_cut)
		return true;
	else
		return false;
}