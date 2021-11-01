#include "stk.h"
#include "cuts.h"

#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpStkContainer.h"
#include "Dmp/DmpPsdContainer.h"
#include "Dmp/DmpStruct.h"

#include "DmpStkTrack.h"
#include "DmpStkTrackHelper.h"

inline bool check_trigger(const std::shared_ptr<DmpEvtHeader> evt_header) {
	auto trigger_mip1 = evt_header->GeneratedTrigger(1);
	auto trigger_mip2 = evt_header->GeneratedTrigger(2);
	auto trigger_HET = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
	auto trigger_LET = evt_header->GeneratedTrigger(4) && evt_header->EnabledTrigger(4);
	auto trigger_MIP = trigger_mip1 || trigger_mip2;
	auto trigger_general = trigger_MIP || trigger_HET || trigger_LET;
    return trigger_general;
}

inline void link_ladders(std::vector<int> &LadderToLayer) {
	for (int ilad = 0; ilad < nSTKladders; ++ilad) {
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
	const std::vector<double> bgoRec_intercept) {

        // Build bgoRecDirection TVector3
        TVector3 vec_s0_a(bgoRec_intercept[0], bgoRec_intercept[1], 0.);
        TVector3 vec_s1_a(bgoRec_intercept[0] + bgoRec_slope[0], bgoRec_intercept[1] + bgoRec_slope[1], 1.);
        bgoRecDirection = (vec_s1_a - vec_s0_a).Unit(); //uni vector pointing from front to back

        // Build bgoRecEntrance TVector3
        double topZ = BGO_TopZ;
        double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
        double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

        if (fabs(topX) > BGO_SideXY || fabs(topY) > BGO_SideXY) {
            // possibly enter from the x-sides
            if (fabs(topX) > BGO_SideXY) {
                if (topX > 0)
                    topX = BGO_SideXY;
                else
                    topX = -BGO_SideXY;
                topZ = (topX - bgoRec_intercept[0]) / bgoRec_slope[0];
                topY = bgoRec_slope[1] * topZ + bgoRec_intercept[1];
                // possibly enter from the y-sides
                if (fabs(topY) > BGO_SideXY) {
                    if (topY > 0)
                        topY = BGO_SideXY;
                    else
                        topY = -BGO_SideXY;
                    topZ = (topY - bgoRec_intercept[1]) / bgoRec_slope[1];
                    topX = bgoRec_slope[0] * topZ + bgoRec_intercept[0];
                }
            }
            //enter from the y-sides
            else if (fabs(topY) > BGO_SideXY) {
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
	const bool best_track = false) {

        std::vector<int> prevHole(2, -2);
        std::vector<int> firstLayer(2, -1);
        std::vector<int> lastLayer(2, -1);
        std::vector<int> lastPoint(2, -1);
        std::vector<unsigned int> track_nHoles_cont(2, 0);

        // Loop on track points to find last layer values
        for (int ip = track->GetNPoints() - 1; ip >= 0; --ip) {
            if (lastLayer[0] == -1) {
                if (track->getHitMeasX(ip) > -99999) {
                    lastPoint[0] = ip;
                    DmpStkSiCluster *cluster = track->GetClusterX(ip, stkclusters.get());
                    auto hardID = cluster->getLadderHardware();
                    lastLayer[0] = LadderToLayer[hardID];
                }
            }
            if (lastLayer[1] == -1) {
                if (track->getHitMeasY(ip) > -99999) {
                    lastPoint[1] = ip;
                    DmpStkSiCluster *cluster = track->GetClusterY(ip, stkclusters.get());
                    auto hardID = cluster->getLadderHardware();
                    lastLayer[1] = LadderToLayer[hardID];
                }
            }
        }

        // Found the number of holes on both X and Y
        for (int ip = 0; ip <= lastPoint[0]; ++ip) {
            if (track->getHitMeasX(ip) > -99999) {
                DmpStkSiCluster *cluster = track->GetClusterX(ip, stkclusters.get());
                auto hardID = cluster->getLadderHardware();
                if (firstLayer[0] == -1)
                    firstLayer[0] = LadderToLayer[hardID];
            }
            else {
                if (firstLayer[0] != -1)
                    ++track_nHoles[0];
                if (ip == prevHole[0] + 1)
                    ++track_nHoles_cont[0];
                prevHole[0] = ip;
            }
        }

        for (int ip = 0; ip <= lastPoint[1]; ++ip) {
            if (track->getHitMeasY(ip) > -99999) {
                DmpStkSiCluster *cluster = track->GetClusterY(ip, stkclusters.get());
                auto hardID = cluster->getLadderHardware();
                if (firstLayer[1] == -1)
                    firstLayer[1] = LadderToLayer[hardID];
            }
            else {
                if (firstLayer[1] != -1)
                    ++track_nHoles[1];
                if (ip == prevHole[1] + 1)
                    ++track_nHoles_cont[1];
                prevHole[1] = ip;
            }
        }

        if (best_track) {
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

inline bool track_selection(
	const std::shared_ptr<DmpEvtBgoRec> bgorec,
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtBgoHits> bgohits,
	const std::shared_ptr<TClonesArray> stkclusters,
	const std::shared_ptr<TClonesArray> stktracks,
    best_track &event_best_track,
    const double evt_corr_energy_gev,
    std::shared_ptr<histos> ps_histos,
    const bool lastcut = false) {
        
        TVector3 bgoRecEntrance;
        TVector3 bgoRecDirection;
        std::vector<int> LadderToLayer(nSTKladders, -1);
        std::vector<DmpStkTrack *> selectedTracks;

        auto weight {ps_histos->GetWeight()};

        link_ladders(LadderToLayer);
        fill_BGO_vectors(
            bgoRecEntrance,
            bgoRecDirection,
            bgoRec_slope,
            bgoRec_intercept);

        // Loop on the tracks
        for (int trIdx = 0; trIdx < stktracks->GetLast() + 1; ++trIdx) {
            std::vector<int> track_nHoles(2, 0);
            std::vector<double> track_slope(2, 0);
            std::vector<double> track_intercept(2, 0);
            std::vector<double> extr_BGO_top(2, 0);

            // *********************

            // Get the track
            auto track = static_cast<DmpStkTrack *>(stktracks->ConstructedAt(trIdx));

            // Get the X and Y clusters
            if (!lastcut) {
                ps_histos->h_STK_X_clusters->Fill(track->getNhitX(), weight);
                ps_histos->h_STK_Y_clusters->Fill(track->getNhitY(), weight);
                ps_histos->h_STK_X_clusters_vs_energy->Fill(evt_corr_energy_gev, track->getNhitX(), weight);
                ps_histos->h_STK_Y_clusters_vs_energy->Fill(evt_corr_energy_gev, track->getNhitY(), weight);
            }
            else {
                ps_histos->h_STK_X_clusters_lastcut->Fill(track->getNhitX(), weight);
                ps_histos->h_STK_Y_clusters_lastcut->Fill(track->getNhitY(), weight);
                ps_histos->h_STK_X_clusters_vs_energy_lastcut->Fill(evt_corr_energy_gev, track->getNhitX(), weight);
                ps_histos->h_STK_Y_clusters_vs_energy_lastcut->Fill(evt_corr_energy_gev, track->getNhitY(), weight);
            }

            get_track_points(
                track,
                stkclusters,
                LadderToLayer,
                track_nHoles,
                event_best_track);

            // Get the X and Y holes
            if (!lastcut) {
                ps_histos->h_STK_X_holes->Fill(track_nHoles[0]);
                ps_histos->h_STK_Y_holes->Fill(track_nHoles[1]);
                ps_histos->h_STK_X_holes_vs_energy->Fill(evt_corr_energy_gev, track_nHoles[0], weight);
                ps_histos->h_STK_Y_holes_vs_energy->Fill(evt_corr_energy_gev, track_nHoles[1], weight);
            }
            else {
                ps_histos->h_STK_X_holes_lastcut->Fill(track_nHoles[0]);
                ps_histos->h_STK_Y_holes_lastcut->Fill(track_nHoles[1]);
                ps_histos->h_STK_X_holes_vs_energy_lastcut->Fill(evt_corr_energy_gev, track_nHoles[0], weight);
                ps_histos->h_STK_Y_holes_vs_energy_lastcut->Fill(evt_corr_energy_gev, track_nHoles[1], weight);
            }

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

            if (!lastcut) {
                ps_histos->h_STK_BGO_TOP_spatial_difference->Fill(drTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_X_difference->Fill(dxTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_Y_difference->Fill(dyTop, weight);
                ps_histos->h_STK_BGO_track_angular_difference->Fill(dAngleTrackBgoRec, weight);

                if (track->getNhitX() == 3 && track->getNhitY() == 3) {
                    ps_histos->h_STK_BGO_TOP_spatial_difference_3_clusters->Fill(drTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_X_difference_3_clusters->Fill(dxTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_Y_difference_3_clusters->Fill(dyTop, weight);
                    ps_histos->h_STK_BGO_track_angular_difference_3_clusters->Fill(dAngleTrackBgoRec, weight);
                }
                else if (track->getNhitX() == 4 && track->getNhitY() == 4) {
                    ps_histos->h_STK_BGO_TOP_spatial_difference_4_clusters->Fill(drTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_X_difference_4_clusters->Fill(dxTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_Y_difference_4_clusters->Fill(dyTop, weight);
                    ps_histos->h_STK_BGO_track_angular_difference_4_clusters->Fill(dAngleTrackBgoRec, weight);
                }
                else if (track->getNhitX() == 5 && track->getNhitY() == 5) {
                    ps_histos->h_STK_BGO_TOP_spatial_difference_5_clusters->Fill(drTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_X_difference_5_clusters->Fill(dxTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_Y_difference_5_clusters->Fill(dyTop, weight);
                    ps_histos->h_STK_BGO_track_angular_difference_5_clusters->Fill(dAngleTrackBgoRec, weight);
                }

            }
            else {
                ps_histos->h_STK_BGO_TOP_spatial_difference_lastcut->Fill(drTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_X_difference_lastcut->Fill(dxTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_Y_difference_lastcut->Fill(dyTop, weight);
                ps_histos->h_STK_BGO_track_angular_difference_lastcut->Fill(dAngleTrackBgoRec, weight);

                if (track->getNhitX() == 3 && track->getNhitY() == 3) {
                    ps_histos->h_STK_BGO_TOP_spatial_difference_3_clusters_lastcut->Fill(drTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_X_difference_3_clusters_lastcut->Fill(dxTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_Y_difference_3_clusters_lastcut->Fill(dyTop, weight);
                    ps_histos->h_STK_BGO_track_angular_difference_3_clusters_lastcut->Fill(dAngleTrackBgoRec, weight);
                }
                else if (track->getNhitX() == 4 && track->getNhitY() == 4) {
                    ps_histos->h_STK_BGO_TOP_spatial_difference_4_clusters_lastcut->Fill(drTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_X_difference_4_clusters_lastcut->Fill(dxTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_Y_difference_4_clusters_lastcut->Fill(dyTop, weight);
                    ps_histos->h_STK_BGO_track_angular_difference_4_clusters_lastcut->Fill(dAngleTrackBgoRec, weight);
                }
                else if (track->getNhitX() == 5 && track->getNhitY() == 5) {
                    ps_histos->h_STK_BGO_TOP_spatial_difference_5_clusters_lastcut->Fill(drTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_X_difference_5_clusters_lastcut->Fill(dxTop, weight);
                    ps_histos->h_STK_BGO_TOP_spatial_Y_difference_5_clusters_lastcut->Fill(dyTop, weight);
                    ps_histos->h_STK_BGO_track_angular_difference_5_clusters_lastcut->Fill(dAngleTrackBgoRec, weight);
                }
            }
            
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
        }
    }

inline bool stk_charge(
    const std::shared_ptr<TClonesArray> stkclusters, 
    best_track &event_best_track, 
    std::shared_ptr<histos> ps_histos,
    const bool lastcut = false) {
    
    auto weight {ps_histos->GetWeight()};

    double cluster_chargeX = -999;
	double cluster_chargeY = -999;

	// Charge correction
	auto track_correction = (event_best_track.myBestTrack).getDirection().CosTheta();

	// Compute charges
	for (auto clIdx = 0; clIdx < event_best_track.n_points; ++clIdx) {
		auto cluster_x = (event_best_track.myBestTrack).GetClusterX(clIdx, stkclusters.get());
		auto cluster_y = (event_best_track.myBestTrack).GetClusterY(clIdx, stkclusters.get());
		if (cluster_x && !cluster_x->getPlane())
			cluster_chargeX = sqrt(cluster_x->getEnergy() * track_correction);
		if (cluster_y && !cluster_y->getPlane())
			cluster_chargeY = sqrt(cluster_y->getEnergy() * track_correction);
	}

	// Check charges
	if (cluster_chargeX != -999 && cluster_chargeY != -999) {
        if (!lastcut) {
            ps_histos->h_STK_charge_X->Fill(cluster_chargeX, weight);
            ps_histos->h_STK_charge_Y->Fill(cluster_chargeY, weight);
            ps_histos->h_STK_charge->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
            ps_histos->h_STK_charge_2D->Fill(cluster_chargeX, cluster_chargeY, weight);

            if (event_best_track.myBestTrack.getNhitX() == 3 && event_best_track.myBestTrack.getNhitY() == 3) {
                ps_histos->h_STK_charge_X_3_clusters->Fill(cluster_chargeX, weight);
                ps_histos->h_STK_charge_Y_3_clusters->Fill(cluster_chargeY, weight);
                ps_histos->h_STK_charge_3_clusters->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
                ps_histos->h_STK_charge_2D_3_clusters->Fill(cluster_chargeX, cluster_chargeY, weight);
            }
            else if (event_best_track.myBestTrack.getNhitX() == 4 && event_best_track.myBestTrack.getNhitY() == 4) {
                ps_histos->h_STK_charge_X_4_clusters->Fill(cluster_chargeX, weight);
                ps_histos->h_STK_charge_Y_4_clusters->Fill(cluster_chargeY, weight);
                ps_histos->h_STK_charge_4_clusters->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
                ps_histos->h_STK_charge_2D_4_clusters->Fill(cluster_chargeX, cluster_chargeY, weight);
            }
            else if (event_best_track.myBestTrack.getNhitX() == 5 && event_best_track.myBestTrack.getNhitY() == 5) {
                ps_histos->h_STK_charge_X_5_clusters->Fill(cluster_chargeX, weight);
                ps_histos->h_STK_charge_Y_5_clusters->Fill(cluster_chargeY, weight);
                ps_histos->h_STK_charge_5_clusters->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
                ps_histos->h_STK_charge_2D_5_clusters->Fill(cluster_chargeX, cluster_chargeY, weight);
            }
        }
        else {
            ps_histos->h_STK_charge_X_lastcut->Fill(cluster_chargeX, weight);
            ps_histos->h_STK_charge_Y_lastcut->Fill(cluster_chargeY, weight);
            ps_histos->h_STK_charge_lastcut->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
            ps_histos->h_STK_charge_2D_lastcut->Fill(cluster_chargeX, cluster_chargeY, weight);

            if (event_best_track.myBestTrack.getNhitX() == 3 && event_best_track.myBestTrack.getNhitY() == 3) {
                ps_histos->h_STK_charge_X_3_clusters_lastcut->Fill(cluster_chargeX, weight);
                ps_histos->h_STK_charge_Y_3_clusters_lastcut->Fill(cluster_chargeY, weight);
                ps_histos->h_STK_charge_3_clusters_lastcut->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
                ps_histos->h_STK_charge_2D_3_clusters_lastcut->Fill(cluster_chargeX, cluster_chargeY, weight);
            }
            else if (event_best_track.myBestTrack.getNhitX() == 4 && event_best_track.myBestTrack.getNhitY() == 4) {
                ps_histos->h_STK_charge_X_4_clusters_lastcut->Fill(cluster_chargeX, weight);
                ps_histos->h_STK_charge_Y_4_clusters_lastcut->Fill(cluster_chargeY, weight);
                ps_histos->h_STK_charge_4_clusters_lastcut->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
                ps_histos->h_STK_charge_2D_4_clusters_lastcut->Fill(cluster_chargeX, cluster_chargeY, weight);
            }
            else if (event_best_track.myBestTrack.getNhitX() == 5 && event_best_track.myBestTrack.getNhitY() == 5) {
                ps_histos->h_STK_charge_X_5_clusters_lastcut->Fill(cluster_chargeX, weight);
                ps_histos->h_STK_charge_Y_5_clusters_lastcut->Fill(cluster_chargeY, weight);
                ps_histos->h_STK_charge_5_clusters_lastcut->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
                ps_histos->h_STK_charge_2D_5_clusters_lastcut->Fill(cluster_chargeX, cluster_chargeY, weight);
            }
        }

    }
}

void stk_distributions(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header, 
    std::shared_ptr<TClonesArray> stkclusters, 
    std::shared_ptr<TClonesArray> stktracks,
    const double evt_energy, 
    const double evt_corr_energy,
    const double evt_energy_gev, 
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos) {

        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
        std::unique_ptr<DmpStkContainer> stkVault = std::make_unique<DmpStkContainer>();
        
        best_track event_best_track;

        double bgo_layer_min_energy     {0};    // Minimum energy per BGO layer
        double bgo_max_energy_ratio     {0.35}; // Maximum energy ratio per layer
        double bgo_shower_axis_delta    {280};  // BGO maximum shower axis delta (mm)
        double bgo_shower_width         {100};  // BGO maximum shower width (mm)
        double STK_BGO_delta_position   {40};   // Linear distance between STK and BGO projections
        double STK_BGO_delta_track      {10};   // Angular distance between BGO/STK tracks (deg)
        int track_X_clusters            {4};    // Number of requested clusters on X tracks
        int track_Y_clusters            {4};    // Number of requested clusters on Y tracks

        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), bgo_layer_min_energy);
        stkVault->scanSTKHits(stkclusters);

        if (check_trigger(evt_header)) {

            auto maxelayer_cut = maxElayer_cut(bgoVault->GetLayerEnergies(), bgo_max_energy_ratio, evt_energy);
            auto maxbarlayer_cut = maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
            auto bgotrack_cut = BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), bgo_shower_axis_delta);

            auto bgofiducial_cut = maxelayer_cut && maxbarlayer_cut && bgotrack_cut;

            if (bgofiducial_cut) {
                auto nbarlayer13_cut = nBarLayer13_cut(bgohits, bgoVault->GetSingleLayerBarNumber(13), evt_energy);
                if (nbarlayer13_cut) {
                    auto maxrms_cut = maxRms_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetRmsLayer(), evt_energy, bgo_shower_width);

                    if (maxrms_cut) {
                        track_selection(
                            bgorec, bgoVault->GetBGOslope(), 
                            bgoVault->GetBGOintercept(), 
                            bgohits, 
                            stkclusters, 
                            stktracks, 
                            event_best_track, 
                            evt_corr_energy_gev, 
                            ps_histos);

                        auto trackselection_cut = track_selection_cut(
                            bgorec, 
                            bgoVault->GetBGOslope(), 
                            bgoVault->GetBGOintercept(), 
                            bgohits, 
                            stkclusters, 
                            stktracks, 
                            event_best_track,
                            STK_BGO_delta_position,
                            STK_BGO_delta_track,
                            track_X_clusters,
                            track_Y_clusters);

                            if (trackselection_cut) {

                                auto weight {ps_histos->GetWeight()};

                                ps_histos->h_BGOrec_sumRms_flast_after_track_selection->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                if (evt_corr_energy_gev>=20 && evt_corr_energy_gev<100)
                                    ps_histos->h_BGOrec_sumRms_flast_after_track_selection_20_100->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=100 && evt_corr_energy_gev<250)
                                    ps_histos->h_BGOrec_sumRms_flast_after_track_selection_100_250->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=250 && evt_corr_energy_gev<500)
                                    ps_histos->h_BGOrec_sumRms_flast_after_track_selection_250_500->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=500 && evt_corr_energy_gev<1000)
                                    ps_histos->h_BGOrec_sumRms_flast_after_track_selection_500_1000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=1000 && evt_corr_energy_gev<3000)
                                    ps_histos->h_BGOrec_sumRms_flast_after_track_selection_1000_3000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=3000 && evt_corr_energy_gev<5000)
                                    ps_histos->h_BGOrec_sumRms_flast_after_track_selection_3000_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=5000)
                                    ps_histos->h_BGOrec_sumRms_flast_after_track_selection_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);

                                stk_charge(stkclusters, event_best_track, ps_histos);
                            }
                    }
                }
            }
        }
    }

void stk_distributions_lastcut(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header, 
    std::shared_ptr<TClonesArray> stkclusters, 
    std::shared_ptr<TClonesArray> stktracks,
    std::shared_ptr<DmpEvtPsdHits> psdhits,
    const double evt_energy, 
    const double evt_corr_energy,
    const double evt_energy_gev, 
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos) {

        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
        std::unique_ptr<DmpStkContainer> stkVault = std::make_unique<DmpStkContainer>();
        std::unique_ptr<DmpPsdContainer> psdVault = std::make_unique<DmpPsdContainer>();

        best_track event_best_track;
        psd_cluster_match clu_matching;

        double bgo_layer_min_energy     {0};    // Minimum energy per BGO layer
        double bgo_max_energy_ratio     {0.35}; // Maximum energy ratio per layer
        double bgo_shower_axis_delta    {280};  // BGO maximum shower axis delta (mm)
        double bgo_shower_width         {100};  // BGO maximum shower width (mm)
        double STK_BGO_delta_position   {40};   // Linear distance between STK and BGO projections
        double STK_BGO_delta_track      {10};   // Angular distance between BGO/STK tracks (deg)
        int track_X_clusters            {4};    // Number of requested clusters on X tracks
        int track_Y_clusters            {4};    // Number of requested clusters on Y tracks
        double psd_min_energy           {0};    // Minimum energy per PSD bar
        double PSD_sharge_sum           {10};   // Sum of PSD charges on X and Y views
        double PSD_single_charge        {2.6};  // PSD charge cut on single view
        
        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), bgo_layer_min_energy);
        stkVault->scanSTKHits(stkclusters);
        psdVault->scanPSDHits(psdhits, psd_min_energy);

        if (check_trigger(evt_header)) {

            auto maxelayer_cut = maxElayer_cut(bgoVault->GetLayerEnergies(), bgo_max_energy_ratio, evt_energy);
            auto maxbarlayer_cut = maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
            auto bgotrack_cut = BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), bgo_shower_axis_delta);

            auto bgofiducial_cut = maxelayer_cut && maxbarlayer_cut && bgotrack_cut;

            if (bgofiducial_cut) {
                auto nbarlayer13_cut = nBarLayer13_cut(bgohits, bgoVault->GetSingleLayerBarNumber(13), evt_energy);
                if (nbarlayer13_cut) {
                    auto maxrms_cut = maxRms_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetRmsLayer(), evt_energy, bgo_shower_width);

                    if (maxrms_cut) {       
                        auto psdcharge_cut = psd_charge_cut(
                            psdVault->getPsdClusterMaxE(),
                            event_best_track,
                            clu_matching,
                            PSD_sharge_sum,
                            PSD_single_charge);

                        if (psdcharge_cut) {

                            track_selection(
                                bgorec, bgoVault->GetBGOslope(), 
                                bgoVault->GetBGOintercept(), 
                                bgohits, 
                                stkclusters, 
                                stktracks, 
                                event_best_track, 
                                evt_corr_energy_gev, 
                                ps_histos,
                                true);
                            
                            // Extract the best track for the event
                            auto trackselection_cut = track_selection_cut(
                                bgorec, 
                                bgoVault->GetBGOslope(), 
                                bgoVault->GetBGOintercept(), 
                                bgohits, 
                                stkclusters, 
                                stktracks, 
                                event_best_track,
                                STK_BGO_delta_position,
                                STK_BGO_delta_track,
                                track_X_clusters,
                                track_Y_clusters);
                            
                            if (trackselection_cut) {
                                stk_charge(stkclusters, event_best_track, ps_histos, true);
                            }
                        }
                    }
                }
            }
        }
    }