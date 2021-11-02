#include "DmpFilterContainer.h"
#include "xtrX_computation.h"

void DmpFilterContainer::Pipeline(
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
	// **** BGO Fiducial Volume ****
	// maxElayer_cut
	output.BGO_fiducial_maxElayer_cut =
		maxElayer_cut(
			bgoVault.GetLayerEnergies(),
			cuts,
			bgoTotalE);
	// maxBarLayer_cut
	output.BGO_fiducial_maxBarLayer_cut =
		maxBarLayer_cut(
			bgoVault.GetLayerBarNumber(),
			bgoVault.GetiMaxLayer(),
			bgoVault.GetIdxBarMaxLayer());
	// BGOTrackContainment_cut
	output.BGO_fiducial_BGOTrackContainment_cut =
		BGOTrackContainment_cut(
			bgoVault.GetBGOslope(),
			bgoVault.GetBGOintercept(),
			cuts);

	output.BGO_fiducial = output.BGO_fiducial_maxElayer_cut && output.BGO_fiducial_maxBarLayer_cut && output.BGO_fiducial_BGOTrackContainment_cut;

	if (output.BGO_fiducial)
	{
		output.all_cut = output.BGO_fiducial;
		// **** nBarLayer13 cut ****
		if (acuts.nBarLayer13)
		{
			output.nBarLayer13_cut = nBarLayer13_cut(
				bgohits,
				bgoVault.GetSingleLayerBarNumber(13),
				bgoTotalE);
			output.all_cut *= output.nBarLayer13_cut;
		}

		if (output.all_cut)
		{
			// **** maxRms cut ****
			if (acuts.maxRms)
			{
				output.maxRms_cut = maxRms_cut(
					bgoVault.GetELayer(),
					bgoVault.GetRmsLayer(),
					bgoTotalE,
					cuts);
				output.all_cut *= output.maxRms_cut;
			}

			if (output.all_cut)
			{
				// **** track selection cut ****
				if (acuts.track_selection)
				{
					output.track_selection_cut = track_selection_cut(
						bgorec,
						bgoVault.GetBGOslope(),
						bgoVault.GetBGOintercept(),
						bgohits,
						stkclusters,
						stktracks,
						cuts);
					output.all_cut *= output.track_selection_cut;
				}

				if (output.all_cut)
				{
					// **** psd stk match cut ****
					if (acuts.psd_stk_match)
					{
						output.psd_stk_match_cut = psd_stk_match_cut(
							bgoVault.GetBGOslope(),
							bgoVault.GetBGOintercept(),
							cuts,
							psdVault.getPsdClusterIdxBegin(),
							psdVault.getPsdClusterZ(),
							psdVault.getPsdClusterMaxECoo());
						output.all_cut *= output.psd_stk_match_cut;
					}

					if (output.all_cut)
					{
						// **** psd charge cut ****
						if (acuts.psd_charge)
						{
							output.psd_charge_measurement = true;
							output.psd_charge_cut = psd_charge_cut(
								psdVault.getPsdClusterMaxE(),
								psdVault.getPsdClusterIdxMaxE(),
								psdVault.getHitZ(),
								psdVault.getGlobalBarID(),
								cuts);
							output.all_cut *= output.psd_charge_cut;
						}

						// **** stk charge cut ****
						if (acuts.stk_charge)
						{
							output.stk_charge_measurement = true;
							output.stk_charge_cut = stk_charge_cut(
								stkclusters,
								cuts);
							output.all_cut *= output.stk_charge_cut;
						}

						if (output.all_cut)
						{
							classifier.xtr = xtrX_computation(
								bgoVault.GetSumRMS(),
								bgoVault.GetSingleFracLayer(13));
							classifier.xtrl = xtrX_computation(
								bgoVault.GetSumRMS(),
								bgoVault.GetSingleFracLayer(bgoVault.GetLastEnergyLayer()));

							output.xtrl_tight_cut = xtrl_tight_cut(classifier.xtrl);
							output.xtrl_loose_cut = xtrl_loose_cut(classifier.xtrl, bgoTotalE_corr);
							++particle_counter.selected_events;
						}
					}
				}
			}
		}
	}
}

const bool DmpFilterContainer::geometric_cut(const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries)
{
	bool passed_geometric_cut = false;

	TVector3 orgPosition;
	orgPosition.SetX(simu_primaries->pv_x);
	orgPosition.SetY(simu_primaries->pv_y);
	orgPosition.SetZ(simu_primaries->pv_z);

#if 0
	// **** Directions Cosines Method

	TVector3 dCos;
	dCos.SetX(simu_primaries->pvpart_cosx);
	dCos.SetY(simu_primaries->pvpart_cosy);
	dCos.SetZ(simu_primaries->pvpart_cosz);

	if (dCos.Z())
	{
		double ratioZ = (BGO_TopZ - orgPosition.Z()) / dCos.Z();
		double actual_X = ratioZ * dCos.X() + orgPosition.X();
		double actual_Y = ratioZ * dCos.Y() + orgPosition.Y();
		if (fabs(actual_X) < BGO_SideXY && fabs(actual_Y) < BGO_SideXY)
			passed_geometric_cut = true;
	}

#else
	// **** Moments Method

	TVector3 orgMomentum;
	orgMomentum.SetX(simu_primaries->pvpart_px);
	orgMomentum.SetY(simu_primaries->pvpart_py);
	orgMomentum.SetZ(simu_primaries->pvpart_pz);

	//auto orgMomentum_theta = orgMomentum.Theta() * TMath::RadToDeg();
	//auto orgMomentum_costheta = cos(orgMomentum.Theta());

	std::vector<double> slope(2, 0);
	std::vector<double> intercept(2, 0);

	slope[0] = orgMomentum.Z() ? orgMomentum.X() / orgMomentum.Z() : -999;
	slope[1] = orgMomentum.Z() ? orgMomentum.Y() / orgMomentum.Z() : -999;
	intercept[0] = orgPosition.X() - slope[0] * orgPosition.Z();
	intercept[1] = orgPosition.Y() - slope[1] * orgPosition.Z();

	double actual_topX = slope[0] * BGO_TopZ + intercept[0];
	double actual_topY = slope[1] * BGO_TopZ + intercept[1];

	double actual_bottomX = slope[0] * BGO_BottomZ + intercept[0];
	double actual_bottomY = slope[1] * BGO_BottomZ + intercept[1];

	if (fabs(actual_topX) < BGO_SideXY && fabs(actual_topY) < BGO_SideXY &&
		fabs(actual_bottomX) < BGO_SideXY && fabs(actual_bottomY) < BGO_SideXY)
		passed_geometric_cut = true;

#endif

	return passed_geometric_cut;
}

const bool DmpFilterContainer::geometric_cut_data(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept)
{
	bool passed_geometric_cut = false;

	double actual_topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
	double actual_topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

	double actual_bottomX = bgoRec_slope[0] * BGO_BottomZ + bgoRec_intercept[0];
	double actual_bottomY = bgoRec_slope[1] * BGO_BottomZ + bgoRec_intercept[1];

	if (fabs(actual_topX) < BGO_SideXY && fabs(actual_topY) < BGO_SideXY &&
		fabs(actual_bottomX) < BGO_SideXY && fabs(actual_bottomY) < BGO_SideXY)
		passed_geometric_cut = true;

	return passed_geometric_cut;
}

const bool DmpFilterContainer::maxElayer_cut(
	const std::vector<double> layer_energies,
	const cuts_conf cuts,
	const double bgoTotalE)
{
	bool passed_maxELayerTotalE_cut = true;

	for (int idxLy = 0; idxLy < DAMPE_bgo_nLayers; ++idxLy)
	{
		auto tmp_ratio = layer_energies[idxLy] / bgoTotalE;
		if (tmp_ratio > cuts.energy_lRatio)
		{
			passed_maxELayerTotalE_cut = false;
			break;
		}
	}

	return passed_maxELayerTotalE_cut;
}

const bool DmpFilterContainer::maxBarLayer_cut(
	const std::vector<std::vector<short>> layerBarNumber,
	const std::vector<int> iMaxLayer,
	const std::vector<int> idxBarMaxLayer)
{
	bool passed_maxBarLayer_cut = true;

	for (auto lIdx = 1; lIdx <= 3; ++lIdx)
	{
		if (layerBarNumber[lIdx].size() == 0)
		{
			passed_maxBarLayer_cut = false;
			break;
		}
		if (iMaxLayer[lIdx] > -1)
			if (idxBarMaxLayer[lIdx] == 0 || idxBarMaxLayer[lIdx] == 21)
			{
				passed_maxBarLayer_cut = false;
				break;
			}
	}

	return passed_maxBarLayer_cut;
}

const bool DmpFilterContainer::BGOTrackContainment_cut(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const cuts_conf cuts)
{
	bool passed_bgo_containment_cut = false;

	double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
	double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

	double bottomX = bgoRec_slope[0] * BGO_BottomZ + bgoRec_intercept[0];
	double bottomY = bgoRec_slope[1] * BGO_BottomZ + bgoRec_intercept[1];

	if (
		fabs(topX) < cuts.shower_axis_delta &&
		fabs(topY) < cuts.shower_axis_delta &&
		fabs(bottomX) < cuts.shower_axis_delta &&
		fabs(bottomY) < cuts.shower_axis_delta)
		passed_bgo_containment_cut = true;

	return passed_bgo_containment_cut;
}

const bool DmpFilterContainer::nBarLayer13_cut(
	const std::shared_ptr<DmpEvtBgoHits> bgohits,
	const std::vector<short> layerBarNumber,
	const double bgoTotalE)
{
	bool passed_nBarLayer13_cut = false;
	unsigned int nTriggeredBGO_13_bars = 0;
	double _GeV = 0.001;

	for (auto it = layerBarNumber.begin(); it != layerBarNumber.end(); ++it)
	{
		auto ihit = *it;
		auto hitE = (bgohits->fEnergy)[ihit];
		if (hitE > 10)
			++nTriggeredBGO_13_bars;
	}
	double nBar13_threshold = 8 * log10(bgoTotalE * _GeV) - 5;
	if (nTriggeredBGO_13_bars < nBar13_threshold)
		passed_nBarLayer13_cut = true;

	return passed_nBarLayer13_cut;
}

const bool DmpFilterContainer::maxRms_cut(
	const std::vector<double> layer_energy,
	const std::vector<double> rmsLayer,
	const double bgoTotalE,
	const cuts_conf cuts)
{
	bool passed_maxRms_cut = false;
	auto max_rms = rmsLayer[0];
	auto eCut = bgoTotalE / 100.;

	for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
		if (layer_energy[lIdx] > eCut)
			if (rmsLayer[lIdx] > max_rms)
				max_rms = rmsLayer[lIdx];

	if (max_rms < cuts.max_rms_shower_width)
		passed_maxRms_cut = true;

	return passed_maxRms_cut;
}

inline void link_ladders(std::vector<int> &LadderToLayer)
{
	for (int ilad = 0; ilad < nSTKladders; ++ilad)
	{
		int iTRB = ilad / 24;
		int iladTRB = ilad % 24;
		int iPlane = 5 - iladTRB / 4;
		int isY = (iTRB / 2 + 1) % 2;
		int iLay = iPlane * 2 + isY;
		LadderToLayer[ilad] = iLay;
	}
}

inline void fill_BGO_vectors(
	TVector3 &bgoRecEntrance,
	TVector3 &bgoRecDirection,
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept)
{
	// Build bgoRecDirection TVector3
	TVector3 vec_s0_a(bgoRec_intercept[0], bgoRec_intercept[1], 0.);
	TVector3 vec_s1_a(bgoRec_intercept[0] + bgoRec_slope[0], bgoRec_intercept[1] + bgoRec_slope[1], 1.);
	bgoRecDirection = (vec_s1_a - vec_s0_a).Unit(); //uni vector pointing from front to back

	// Build bgoRecEntrance TVector3
	double topZ = BGO_TopZ;
	double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
	double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

	if (fabs(topX) > BGO_SideXY || fabs(topY) > BGO_SideXY)
	{
		// possibly enter from the x-sides
		if (fabs(topX) > BGO_SideXY)
		{
			if (topX > 0)
				topX = BGO_SideXY;
			else
				topX = -BGO_SideXY;
			topZ = (topX - bgoRec_intercept[0]) / bgoRec_slope[0];
			topY = bgoRec_slope[1] * topZ + bgoRec_intercept[1];
			// possibly enter from the y-sides
			if (fabs(topY) > BGO_SideXY)
			{
				if (topY > 0)
					topY = BGO_SideXY;
				else
					topY = -BGO_SideXY;
				topZ = (topY - bgoRec_intercept[1]) / bgoRec_slope[1];
				topX = bgoRec_slope[0] * topZ + bgoRec_intercept[0];
			}
		}
		//enter from the y-sides
		else if (fabs(topY) > BGO_SideXY)
		{
			if (topY > 0)
				topY = BGO_SideXY;
			else
				topY = -BGO_SideXY;
			topZ = (topY - bgoRec_intercept[1]) / bgoRec_slope[1];
			topX = bgoRec_slope[0] * topZ + bgoRec_intercept[0];
		}
	}

	bgoRecEntrance[0] = topX;
	bgoRecEntrance[1] = topY;
	bgoRecEntrance[2] = topZ;
}

inline void get_track_points(
	DmpStkTrack *track,
	const std::shared_ptr<TClonesArray> stkclusters,
	const std::vector<int> LadderToLayer,
	std::vector<int> &track_nHoles,
	best_track &event_best_track,
	const bool best_track = false)
{
	std::vector<int> prevHole(2, -2);
	std::vector<int> firstLayer(2, -1);
	std::vector<int> lastLayer(2, -1);
	std::vector<int> lastPoint(2, -1);
	std::vector<unsigned int> track_nHoles_cont(2, 0);

	// Loop on track points to find last layer values
	for (int ip = track->GetNPoints() - 1; ip >= 0; --ip)
	{
		if (lastLayer[0] == -1)
		{
			if (track->getHitMeasX(ip) > -99999)
			{
				lastPoint[0] = ip;
				DmpStkSiCluster *cluster = track->GetClusterX(ip, stkclusters.get());
				auto hardID = cluster->getLadderHardware();
				lastLayer[0] = LadderToLayer[hardID];
			}
		}
		if (lastLayer[1] == -1)
		{
			if (track->getHitMeasY(ip) > -99999)
			{
				lastPoint[1] = ip;
				DmpStkSiCluster *cluster = track->GetClusterY(ip, stkclusters.get());
				auto hardID = cluster->getLadderHardware();
				lastLayer[1] = LadderToLayer[hardID];
			}
		}
	}

	// Found the number of holes on both X and Y
	for (int ip = 0; ip <= lastPoint[0]; ++ip)
	{
		if (track->getHitMeasX(ip) > -99999)
		{
			DmpStkSiCluster *cluster = track->GetClusterX(ip, stkclusters.get());
			auto hardID = cluster->getLadderHardware();
			if (firstLayer[0] == -1)
				firstLayer[0] = LadderToLayer[hardID];
		}
		else
		{
			if (firstLayer[0] != -1)
				++track_nHoles[0];
			if (ip == prevHole[0] + 1)
				++track_nHoles_cont[0];
			prevHole[0] = ip;
		}
	}

	for (int ip = 0; ip <= lastPoint[1]; ++ip)
	{
		if (track->getHitMeasY(ip) > -99999)
		{
			DmpStkSiCluster *cluster = track->GetClusterY(ip, stkclusters.get());
			auto hardID = cluster->getLadderHardware();
			if (firstLayer[1] == -1)
				firstLayer[1] = LadderToLayer[hardID];
		}
		else
		{
			if (firstLayer[1] != -1)
				++track_nHoles[1];
			if (ip == prevHole[1] + 1)
				++track_nHoles_cont[1];
			prevHole[1] = ip;
		}
	}

	if (best_track)
	{
		event_best_track.n_points = track->GetNPoints();
		for (int idx = 0; idx < 2; ++idx)
			event_best_track.n_holes[idx] = track_nHoles[idx];
		event_best_track.track_slope[0] = track->getTrackParams().getSlopeX();
		event_best_track.track_slope[1] = track->getTrackParams().getSlopeY();
		event_best_track.track_intercept[0] = track->getTrackParams().getInterceptX();
		event_best_track.track_intercept[1] = track->getTrackParams().getInterceptY();
		event_best_track.track_direction = (track->getDirection()).Unit();
	}
}

const bool DmpFilterContainer::track_selection_cut(
	const std::shared_ptr<DmpEvtBgoRec> bgorec,
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtBgoHits> bgohits,
	const std::shared_ptr<TClonesArray> stkclusters,
	const std::shared_ptr<TClonesArray> stktracks,
	const cuts_conf cuts)
{
	bool passed_track_selection_cut = false;

	TVector3 bgoRecEntrance;
	TVector3 bgoRecDirection;
	std::vector<int> LadderToLayer(nSTKladders, -1);
	std::vector<DmpStkTrack *> selectedTracks;

	link_ladders(LadderToLayer);
	fill_BGO_vectors(
		bgoRecEntrance,
		bgoRecDirection,
		bgoRec_slope,
		bgoRec_intercept);

	// Loop on the tracks
	for (int trIdx = 0; trIdx < stktracks->GetLast() + 1; ++trIdx)
	{
		std::vector<int> track_nHoles(2, 0);
		std::vector<double> track_slope(2, 0);
		std::vector<double> track_intercept(2, 0);
		std::vector<double> extr_BGO_top(2, 0);

		// *********************

		// Get the track
		auto track = static_cast<DmpStkTrack *>(stktracks->ConstructedAt(trIdx));

		// Reject tracks with not enough X and Y clusters
		if (track->getNhitX() < cuts.track_X_clusters || track->getNhitY() < cuts.track_Y_clusters)
			continue;

		get_track_points(
			track,
			stkclusters,
			LadderToLayer,
			track_nHoles,
			event_best_track);

		if (track_nHoles[0] > 1 || track_nHoles[1] > 1)
			continue;

		// Find slope and intercept
		track_slope[0] = track->getTrackParams().getSlopeX();
		track_slope[1] = track->getTrackParams().getSlopeY();
		track_intercept[0] = track->getTrackParams().getInterceptX();
		track_intercept[1] = track->getTrackParams().getInterceptY();
		TVector3 trackDirection = (track->getDirection()).Unit();

		// Extrapolate to the top of BGO
		for (int coord = 0; coord < 2; ++coord)
			extr_BGO_top[coord] = track_slope[coord] * BGO_TopZ + track_intercept[coord];

		// Evaluate distance between Top STK and BGO points
		double dxTop = extr_BGO_top[0] - bgoRecEntrance[0];
		double dyTop = extr_BGO_top[1] - bgoRecEntrance[1];
		double drTop = sqrt(pow(dxTop, 2) + pow(dyTop, 2));
		// Evaluate angular distance between STK track and BGO Rec track
		double dAngleTrackBgoRec = trackDirection.Angle(bgoRecDirection) * TMath::RadToDeg();

		if (drTop > cuts.STK_BGO_delta_position)
			continue;
		if (dAngleTrackBgoRec > cuts.STK_BGO_delta_track)
			continue;

		selectedTracks.push_back(track);
	}

	// Sort selected tracks vector
	DmpStkTrackHelper tHelper(stktracks.get(), true, bgorec.get(), bgohits.get());
	tHelper.MergeSort(selectedTracks, &DmpStkTrackHelper::TracksCompare);

	if (selectedTracks.size() > 0)
	{
		DmpStkTrack *selected_track = static_cast<DmpStkTrack *>(selectedTracks[0]);
		std::vector<int> track_nHoles(2, 0);

		// Fill best track structure
		get_track_points(
			selected_track,
			stkclusters,
			LadderToLayer,
			track_nHoles,
			event_best_track,
			true);

		event_best_track.extr_BGO_topX = event_best_track.track_slope[0] * BGO_TopZ + event_best_track.track_intercept[0];
		event_best_track.extr_BGO_topY = event_best_track.track_slope[1] * BGO_TopZ + event_best_track.track_intercept[1];

		event_best_track.STK_BGO_topX_distance = event_best_track.extr_BGO_topX - bgoRecEntrance[0];
		event_best_track.STK_BGO_topY_distance = event_best_track.extr_BGO_topY - bgoRecEntrance[1];
		event_best_track.angular_distance_STK_BGO = event_best_track.track_direction.Angle(bgoRecDirection) * TMath::RadToDeg();
		event_best_track.STK_BGO_topY_distance = sqrt(pow(event_best_track.STK_BGO_topX_distance, 2) + pow(event_best_track.STK_BGO_topY_distance, 2));

		event_best_track.myBestTrack = static_cast<DmpStkTrack>(*selected_track);

		passed_track_selection_cut = true;
	}

	return passed_track_selection_cut;
}

const bool DmpFilterContainer::psd_stk_match_cut(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const cuts_conf cuts,
	const std::vector<std::vector<short>> psdCluster_idxBeg,
	const std::vector<std::vector<double>> psdCluster_Z,
	const std::vector<std::vector<double>> psdCluster_maxEcoordinate)
{
	bool passed_stk_psd_match = false;

	for (int nLayer = 0; nLayer < DAMPE_psd_nLayers; ++nLayer)
	{
		for (unsigned int iclu = 0; iclu < psdCluster_idxBeg[nLayer].size(); ++iclu)
		{
			bool IsMeasuringX = nLayer % 2;
			double hitZ = psdCluster_Z[nLayer][iclu];
			double thisCoord = psdCluster_maxEcoordinate[nLayer][iclu];

			// Get distance between actual coordinate and BGO rec coordinate
			double projCoord_bgoRec = IsMeasuringX ? bgoRec_slope[0] * hitZ + bgoRec_intercept[0] : bgoRec_slope[1] * hitZ + bgoRec_intercept[1];
			double dX_bgoRec = thisCoord - projCoord_bgoRec;
			if (fabs(dX_bgoRec) < fabs(clu_matching.dxCloPsdClu_bgoRec[nLayer]))
			{
				clu_matching.dxCloPsdClu_bgoRec[nLayer] = dX_bgoRec;
				clu_matching.icloPsdClu_bgoRec[nLayer] = iclu;
			}

			// Get distance between actual coordinate and best track coordinate
			double projCoord_track = IsMeasuringX ? event_best_track.track_slope[0] * hitZ + event_best_track.track_intercept[0] : event_best_track.track_slope[1] * hitZ + event_best_track.track_intercept[1];
			double dX_track = thisCoord - projCoord_track;

			if (fabs(dX_track) < fabs(clu_matching.dxCloPsdClu_track[nLayer]))
			{
				clu_matching.dxCloPsdClu_track[nLayer] = dX_track;
				clu_matching.icloPsdClu_track[nLayer] = iclu;
			}
		}
	}

	passed_stk_psd_match = (fabs(clu_matching.dxCloPsdClu_track[0]) < cuts.STK_PSD_delta_position && fabs(clu_matching.dxCloPsdClu_track[1]) < cuts.STK_PSD_delta_position) ? true : false;
	return passed_stk_psd_match;
}

const bool DmpFilterContainer::psd_charge_cut(
	const std::vector<std::vector<double>> psdCluster_maxE,
	const std::vector<std::vector<short>> psdCluster_idxMaxE,
	const std::vector<double> hitZ,
	const std::vector<short> globalBarID,
	const cuts_conf cuts)
{
	bool passed_psd_charge_cut = false;
	bool passed_psd_charge_sum = false;
	bool passed_psd_He_cut = false;

	double psd_chargeX, psd_chargeY;
	double psd_chargeX_ecor, psd_chargeY_ecor;

	// Charge correction
	auto track_correction = (event_best_track.myBestTrack).getDirection().CosTheta();

	// Get Y charge
	if (clu_matching.icloPsdClu_track[0] > -1)
	{
		auto energy_ClusterYTrack = psdCluster_maxE[0][clu_matching.icloPsdClu_track[0]];
		auto energy_ClusterYTrack_corr = track_correction * energy_ClusterYTrack;
		psd_chargeY = sqrt(energy_ClusterYTrack_corr);
		auto imax = psdCluster_idxMaxE[0][clu_matching.icloPsdClu_track[0]];
		auto hitZpos = hitZ[imax];
		auto gid = globalBarID[imax];
		auto pos = event_best_track.track_slope[1] * hitZpos + event_best_track.track_intercept[1];
		auto PsdEC_tmp = gPsdECor->GetPsdECor(gid, pos / 10.);
		auto energy_ClusterYPsdECor = energy_ClusterYTrack_corr * PsdEC_tmp;
		psd_chargeY_ecor = sqrt(energy_ClusterYPsdECor * track_correction * PsdEC_tmp);
	}

	// Get X charge
	if (clu_matching.icloPsdClu_track[1] > -1)
	{
		auto energy_ClusterXTrack = psdCluster_maxE[1][clu_matching.icloPsdClu_track[1]];
		auto energy_ClusterXTrack_corr = track_correction * energy_ClusterXTrack;
		psd_chargeX = sqrt(energy_ClusterXTrack_corr);
		auto imax = psdCluster_idxMaxE[1][clu_matching.icloPsdClu_track[1]];
		auto hitZpos = hitZ[imax];
		auto gid = globalBarID[imax];
		auto pos = event_best_track.track_slope[0] * hitZpos + event_best_track.track_intercept[0];
		auto PsdEC_tmp = gPsdECor->GetPsdECor(gid, pos / 10.);
		auto energy_ClusterXPsdECor = energy_ClusterXTrack_corr * PsdEC_tmp;
		psd_chargeX_ecor = sqrt(energy_ClusterXPsdECor * track_correction * PsdEC_tmp);
	}

	// Fill PSD ccharge struct
	extracted_psd_charge.chargeX = psd_chargeX;
	extracted_psd_charge.chargeY = psd_chargeY;

	if ((psd_chargeX + psd_chargeY) < cuts.PSD_charge_sum)
		passed_psd_charge_sum = true;

	if (psd_chargeX < cuts.PSD_charge || psd_chargeY < cuts.PSD_charge)
		passed_psd_He_cut = true;

	passed_psd_charge_cut = passed_psd_charge_sum && passed_psd_He_cut;

	return passed_psd_charge_cut;
}

const bool DmpFilterContainer::stk_charge_cut(
	const std::shared_ptr<TClonesArray> stkclusters,
	const cuts_conf cuts)
{
	bool passed_stk_charge_cut = false;

	double cluster_chargeX = -999;
	double cluster_chargeY = -999;

	// Charge correction
	auto track_correction = (event_best_track.myBestTrack).getDirection().CosTheta();

	// Compute charges
	for (auto clIdx = 0; clIdx < event_best_track.n_points; ++clIdx)
	{
		auto cluster_x = (event_best_track.myBestTrack).GetClusterX(clIdx, stkclusters.get());
		auto cluster_y = (event_best_track.myBestTrack).GetClusterY(clIdx, stkclusters.get());
		if (cluster_x && !cluster_x->getPlane())
			cluster_chargeX = sqrt(cluster_x->getEnergy() * track_correction);
		if (cluster_y && !cluster_y->getPlane())
			cluster_chargeY = sqrt(cluster_y->getEnergy() * track_correction);
	}

	// Check charges
	if (cluster_chargeX == -999 || cluster_chargeY == -999)
		return passed_stk_charge_cut;

	// Fill STK charge struct
	extracted_stk_charge.chargeX = cluster_chargeX;
	extracted_stk_charge.chargeY = cluster_chargeY;

	// Compute mean charge
	auto mean_charge = 0.5 * (cluster_chargeX + cluster_chargeY);

	// Check STK charge to select electrons
	if (mean_charge < cuts.STK_charge)
		passed_stk_charge_cut = true;

	return passed_stk_charge_cut;
}

const bool DmpFilterContainer::xtrl_tight_cut(
	const double input_xtrl,
	const double cut_value)
{
	bool passed_xtrl = false;
	if (input_xtrl < cut_value)
		passed_xtrl = true;
	return passed_xtrl;
}

const bool DmpFilterContainer::xtrl_loose_cut(
	const double input_xtrl,
	const double energy)
{
	bool passed_xtrl = false;
	double _GeV = 0.001;
	double cut_value = 3 * pow(log10(energy * _GeV / 200), 2) + 12;
	if (input_xtrl < cut_value)
		passed_xtrl = true;
	return passed_xtrl;
}

void DmpFilterContainer::Reset()
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

void DmpFilterContainer::reset_stk_best_track()
{
	event_best_track.n_points = 0;
	event_best_track.n_holes = {0, 0};
	event_best_track.track_slope = {-999, -999};
	event_best_track.track_intercept = {-999, -999};
	event_best_track.track_direction.SetXYZ(-999, -999, -999);
	event_best_track.extr_BGO_topX = -999;
	event_best_track.extr_BGO_topY = -999;
	event_best_track.STK_BGO_topX_distance = -999;
	event_best_track.STK_BGO_topY_distance = -999;
	event_best_track.angular_distance_STK_BGO = -999;
	event_best_track.myBestTrack = DmpStkTrack();
}

void DmpFilterContainer::reset_psd_clusters()
{
	clu_matching.icloPsdClu = std::vector<int>(DAMPE_psd_nLayers, -999);
	clu_matching.dxCloPsdClu = std::vector<double>(DAMPE_psd_nLayers, -999);
	clu_matching.icloPsdCluMaxHit = std::vector<int>(DAMPE_psd_nLayers, -999);
	clu_matching.dxCloPsdCluMaxHit = std::vector<double>(DAMPE_psd_nLayers, -999);
	clu_matching.icloPsdClu_bgoRec = std::vector<int>(DAMPE_psd_nLayers, -999);
	clu_matching.dxCloPsdClu_bgoRec = std::vector<double>(DAMPE_psd_nLayers, -999);
	clu_matching.icloPsdClu_track = std::vector<int>(DAMPE_psd_nLayers, -999);
	clu_matching.dxCloPsdClu_track = std::vector<double>(DAMPE_psd_nLayers, -999);
	clu_matching.icloPsdClu2_track = std::vector<int>(DAMPE_psd_nLayers, -999);
	clu_matching.dxCloPsdClu2_track = std::vector<double>(DAMPE_psd_nLayers, -999);
}

void DmpFilterContainer::reset_psd_charges()
{
	extracted_psd_charge.chargeX = -999;
	extracted_psd_charge.chargeY = -999;
}

void DmpFilterContainer::reset_stk_charges()
{
	extracted_stk_charge.chargeX = -999;
	extracted_stk_charge.chargeY = -999;
}

void DmpFilterContainer::reset_classifiers()
{
	classifier.xtr = -999;
	classifier.xtrl = -999;
}

void DmpFilterContainer::reset_filter_output()
{
	output.out_energy_range = false;
	output.geometric_before_trigger = false;
	output.evt_in_saa = false;
	output.trigger_check = false;
	output.evt_triggered = false;
	output.correct_bgo_reco = false;
	output.good_event = false;
	output.geometric = false;
	output.BGO_fiducial = false;
	output.BGO_fiducial_maxElayer_cut = false;
	output.BGO_fiducial_maxBarLayer_cut = false;
	output.BGO_fiducial_BGOTrackContainment_cut = false;
	output.nBarLayer13_cut = false;
	output.maxRms_cut = false;
	output.track_selection_cut = false;
	output.psd_stk_match_cut = false;
	output.psd_charge_cut = false;
	output.stk_charge_cut = false;
	output.psd_charge_measurement = false;
	output.stk_charge_measurement = false;
	output.all_cut = false;
}

void DmpFilterContainer::reset_time()
{
	time.second = 0;
}

void DmpFilterContainer::reset_trigger()
{
	evt_trigger.mip1 = false;
	evt_trigger.mip2 = false;
	evt_trigger.MIP = false;
	evt_trigger.HET = false;
	evt_trigger.LET = false;
	evt_trigger.general = false;
}

void DmpFilterContainer::CheckGeometry(
	const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept)
{
	if (!output.evt_in_saa)
		if (!output.out_energy_range)
		{
			if (simu_primaries!=nullptr)
			{
				auto geocut_result = geometric_cut(simu_primaries);
				!output.trigger_check ? output.geometric_before_trigger = geocut_result : output.geometric = geocut_result;
			}
			else
			{
				auto geocut_result = geometric_cut_data(bgoRec_slope, bgoRec_intercept);
				output.geometric = geocut_result;
			}
		}
}

void DmpFilterContainer::checkBGOreco(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries)
{
	bool bgoreco_status = false;
	if ((bgoRec_slope[0] == 0 && bgoRec_intercept[0] == 0) || (bgoRec_slope[1] == 0 && bgoRec_intercept[1] == 0))
	{
		if (simu_primaries)
		{
			TVector3 orgPosition;
			orgPosition.SetX(simu_primaries->pv_x);
			orgPosition.SetY(simu_primaries->pv_y);
			orgPosition.SetZ(simu_primaries->pv_z);

			TVector3 orgMomentum;
			orgMomentum.SetX(simu_primaries->pvpart_px);
			orgMomentum.SetY(simu_primaries->pvpart_py);
			orgMomentum.SetZ(simu_primaries->pvpart_pz);

			std::vector<double> slope(2, 0);
			std::vector<double> intercept(2, 0);

			slope[0] = orgMomentum.Z() ? orgMomentum.X() / orgMomentum.Z() : -999;
			slope[1] = orgMomentum.Z() ? orgMomentum.Y() / orgMomentum.Z() : -999;
			intercept[0] = orgPosition.X() - slope[0] * orgPosition.Z();
			intercept[1] = orgPosition.Y() - slope[1] * orgPosition.Z();

			double actual_X = slope[0] * BGO_TopZ + intercept[0];
			double actual_Y = slope[1] * BGO_TopZ + intercept[1];
			double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
			double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

			int position_sensitivity = 30;

			if (fabs(actual_X - topX) < position_sensitivity && fabs(actual_Y - topY) < position_sensitivity)
				bgoreco_status = true;
		}
	}
	else
		bgoreco_status = true;
	output.correct_bgo_reco = bgoreco_status;
}

void DmpFilterContainer::check_trigger(const std::shared_ptr<DmpEvtHeader> evt_header)
{
	evt_trigger.mip1 = evt_header->GeneratedTrigger(1);
	evt_trigger.mip2 = evt_header->GeneratedTrigger(2);
	evt_trigger.HET = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
	evt_trigger.LET = evt_header->GeneratedTrigger(4) && evt_header->EnabledTrigger(4);
	evt_trigger.MIP = evt_trigger.mip1 || evt_trigger.mip2;
	evt_trigger.general = evt_trigger.MIP || evt_trigger.HET || evt_trigger.LET;
	if (evt_trigger.HET)
		++particle_counter.triggered_events;
	output.evt_triggered = evt_trigger.HET;
	output.trigger_check = true;
}

const bool DmpFilterContainer::CheckIncomingEvent(
	const std::shared_ptr<DmpEvtHeader> evt_header,
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries)
{
	if (!output.evt_in_saa)
	{
		check_trigger(evt_header);
		checkBGOreco(bgoRec_slope, bgoRec_intercept, simu_primaries);
		if (output.evt_triggered && output.correct_bgo_reco && !output.out_energy_range)
			output.good_event = true;
	}
	return output.good_event;
}

void DmpFilterContainer::EnergyCheck(
	const cuts_conf &cuts,
	const double bgoTotalE_corr,
	const double min_energy,
	const double max_energy)
{
	if (!output.evt_in_saa)
	{
		const double _GeV = 0.001;
		if (bgoTotalE_corr * _GeV < cuts.min_event_energy || bgoTotalE_corr * _GeV > cuts.max_event_energy)
		{
			++particle_counter.events_out_range;
			output.out_energy_range = true;
		}
		else
			++particle_counter.events_in_range;
	}
}

void DmpFilterContainer::UpdateEvtCounter()
{
	++particle_counter.event_counter;
}

const filter_output DmpFilterContainer::GetFilterOutput()
{
	return output;
}

const psd_charge DmpFilterContainer::GetPSDCharge()
{
	return extracted_psd_charge;
}

const stk_charge DmpFilterContainer::GetSTKCharge()
{
	return extracted_stk_charge;
}

const best_track DmpFilterContainer::GetBestTrack()
{
	return event_best_track;
}

DmpStkTrack DmpFilterContainer::GetBestTrackObj()
{
	return event_best_track.myBestTrack;
}

const bgo_classifiers DmpFilterContainer::GetClassifiers()
{
	return classifier;
}

const unsigned int DmpFilterContainer::GetStatEvtCounter()
{
	return particle_counter.event_counter;
}

const unsigned int DmpFilterContainer::GetStatEvtInRange()
{
	return particle_counter.events_in_range;
}

const unsigned int DmpFilterContainer::GetStatEvtOutRange()
{
	return particle_counter.events_out_range;
}

const unsigned int DmpFilterContainer::GetStatiEvtTrigger()
{
	return particle_counter.triggered_events;
}

const unsigned int DmpFilterContainer::GetStatEvtSelection()
{
	return particle_counter.selected_events;
}

const unsigned int DmpFilterContainer::GetStatEvtSAA()
{
	return particle_counter.events_in_saa;
}

void DmpFilterContainer::UpdateEvtTime(
	const std::shared_ptr<TChain> dmpch,
	const std::shared_ptr<DmpEvtHeader> evt_header)
{
	dmpch->GetEntry(0);
	time.start_second = evt_header->GetSecond();
	time.start_msecond = evt_header->GetMillisecond();
	dmpch->GetEntry(dmpch->GetEntries() - 1);
	time.end_second = evt_header->GetSecond();
	time.end_msecond = evt_header->GetMillisecond();
}

void DmpFilterContainer::PrintDataInfo(
	const std::shared_ptr<TChain> dmpch,
	const bool mc)
{
	std::cout << "INFO: \tTotal number of events: " << dmpch->GetEntries() << std::endl;
	if (!mc)
	{
		std::cout << "INFO: \tStart Time (second): " << time.start_second << std::endl;
		std::cout << "INFO: \tStart Time (millisecond): " << time.start_msecond << std::endl;
		std::cout << "INFO: \tEnd Time (second): " << time.end_second << std::endl;
		std::cout << "INFO: \tEnd Time (millisecond): " << time.end_msecond << "\n\n";
	}
}

void DmpFilterContainer::SAACheck(
	const std::shared_ptr<DmpEvtHeader> evt_header,
	const std::shared_ptr<DmpFilterOrbit> pFilter)
{
	// Extract data time information
	time.second = evt_header->GetSecond();
	time.msecond = evt_header->GetMillisecond();
	// SAA
	if (pFilter->IsInSAA(time.second))
	{
		++particle_counter.events_in_saa;
		output.evt_in_saa = true;
	}
}

const data_evt_time DmpFilterContainer::GetDataTime()
{
	return time;
}

const trigger_info DmpFilterContainer::GetTrigger()
{
	return evt_trigger;
}
