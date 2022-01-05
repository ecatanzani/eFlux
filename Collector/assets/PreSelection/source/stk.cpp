#include "stk.h"
#include "cuts.h"
#include "charge.h"
#include "trigger.h"

#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpStkContainer.h"
#include "Dmp/DmpStruct.h"

#include "DmpStkTrackHelper.h"

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
    std::shared_ptr<histos> ps_histos,
    std::shared_ptr<config> cuts_config) {

        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
        std::unique_ptr<DmpStkContainer> stkVault = std::make_unique<DmpStkContainer>();
        
        best_track event_best_track;

        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), (cuts_config->GetCutsConfig()).bgo_layer_min_energy);
        stkVault->scanSTKHits(stkclusters);

        if (check_trigger(evt_header)) {

            auto maxelayer_cut = maxElayer_cut(bgoVault->GetLayerEnergies(), (cuts_config->GetCutsConfig()).bgo_max_energy_ratio, evt_energy);
            auto maxbarlayer_cut = maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
            auto bgotrack_cut = BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), (cuts_config->GetCutsConfig()).bgo_shower_axis_delta);

            auto bgofiducial_cut = maxelayer_cut && maxbarlayer_cut && bgotrack_cut;

            if (bgofiducial_cut) {
                auto nbarlayer13_cut = nBarLayer13_cut(bgohits, bgoVault->GetSingleLayerBarNumber(13), evt_energy);
                if (nbarlayer13_cut) {
                    auto maxrms_cut = maxRms_cut(bgoVault->GetELayer(), bgoVault->GetRmsLayer(), evt_energy, (cuts_config->GetCutsConfig()).bgo_shower_width);

                    if (maxrms_cut) {

                        track(
                            bgorec, bgoVault->GetBGOslope(), 
                            bgoVault->GetBGOintercept(), 
                            bgohits, 
                            stkclusters, 
                            stktracks,
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
                            (cuts_config->GetCutsConfig()).STK_BGO_delta_position,
                            (cuts_config->GetCutsConfig()).STK_BGO_delta_track,
                            (cuts_config->GetCutsConfig()).track_X_clusters,
                            (cuts_config->GetCutsConfig()).track_Y_clusters,
                            (cuts_config->GetCutsConfig()).track_X_holes,
                            (cuts_config->GetCutsConfig()).track_Y_holes);

                        auto weight {ps_histos->GetWeight()};

                        if (trackselection_cut) {
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

                            ps_histos->h_trackselection_lastcut_pass->Fill(evt_corr_energy_gev, weight);
                        }
                        else
                            ps_histos->h_trackselection_lastcut_fail->Fill(evt_corr_energy_gev, weight);
                        ps_histos->h_trackselection_lastcut->Fill(evt_corr_energy_gev, weight);
                    }
                }
            }
        }
    }

void BGO_vectors(
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

void ladders(std::vector<int> &LadderToLayer) {
	for (int ilad = 0; ilad < nSTKladders; ++ilad) {
		int iTRB = ilad / 24;
		int iladTRB = ilad % 24;
		int iPlane = 5 - iladTRB / 4;
		int isY = (iTRB / 2 + 1) % 2;
		int iLay = iPlane * 2 + isY;
		LadderToLayer[ilad] = iLay;
	}
}

void track_points(
	DmpStkTrack *track,
	const std::shared_ptr<TClonesArray> stkclusters,
	const std::vector<int> LadderToLayer,
	std::vector<int> &track_nHoles) {

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
    }

void track(
	const std::shared_ptr<DmpEvtBgoRec> bgorec,
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtBgoHits> bgohits,
	const std::shared_ptr<TClonesArray> stkclusters,
	const std::shared_ptr<TClonesArray> stktracks,
    const double evt_corr_energy_gev,
    std::shared_ptr<histos> ps_histos) {
        
        TVector3 bgoRecEntrance;
        TVector3 bgoRecDirection;
        std::vector<int> LadderToLayer(nSTKladders, -1);
        std::vector<DmpStkTrack *> selectedTracks;

        auto weight {ps_histos->GetWeight()};

        ladders(LadderToLayer);
        BGO_vectors(bgoRecEntrance, bgoRecDirection, bgoRec_slope, bgoRec_intercept);

        ps_histos->h_STK_tracks->Fill(stktracks->GetLast() + 1, weight);
        ps_histos->h_STK_tracks_vs_energy->Fill(evt_corr_energy_gev, stktracks->GetLast() + 1, weight);

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
            ps_histos->h_STK_X_clusters->Fill(track->getNhitX(), weight);
            ps_histos->h_STK_Y_clusters->Fill(track->getNhitY(), weight);
            ps_histos->h_STK_X_clusters_vs_energy->Fill(evt_corr_energy_gev, track->getNhitX(), weight);
            ps_histos->h_STK_Y_clusters_vs_energy->Fill(evt_corr_energy_gev, track->getNhitY(), weight);

            track_points(
                track,
                stkclusters,
                LadderToLayer,
                track_nHoles);

            // Get the X and Y holes
            ps_histos->h_STK_X_holes->Fill(track_nHoles[0]);
            ps_histos->h_STK_Y_holes->Fill(track_nHoles[1]);
            ps_histos->h_STK_X_holes_vs_energy->Fill(evt_corr_energy_gev, track_nHoles[0], weight);
            ps_histos->h_STK_Y_holes_vs_energy->Fill(evt_corr_energy_gev, track_nHoles[1], weight);

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

            selectedTracks.push_back(track);

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
            else if (track->getNhitX() == 6 && track->getNhitY() == 6) {
                ps_histos->h_STK_BGO_TOP_spatial_difference_6_clusters->Fill(drTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_X_difference_6_clusters->Fill(dxTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_Y_difference_6_clusters->Fill(dyTop, weight);
                ps_histos->h_STK_BGO_track_angular_difference_6_clusters->Fill(dAngleTrackBgoRec, weight);
            }
        }

        // Sort selected tracks vector
        DmpStkTrackHelper tHelper(stktracks.get(), true, bgorec.get(), bgohits.get());
        tHelper.MergeSort(selectedTracks, &DmpStkTrackHelper::TracksCompare);

        if (selectedTracks.size() > 0) {
            DmpStkTrack *selected_track = static_cast<DmpStkTrack *>(selectedTracks[0]);

            if (selected_track->getNhitX() == 3 && selected_track->getNhitY() == 3)
                ps_histos->h_STK_best_track_clusters->Fill(selected_track->getNhitX(), weight);
            else if (selected_track->getNhitX() == 4 && selected_track->getNhitY() == 4)
                ps_histos->h_STK_best_track_clusters->Fill(selected_track->getNhitX(), weight);
            else if (selected_track->getNhitX() == 5 && selected_track->getNhitY() == 5)
                ps_histos->h_STK_best_track_clusters->Fill(selected_track->getNhitX(), weight);
            else if (selected_track->getNhitX() == 6 && selected_track->getNhitY() == 6)
                ps_histos->h_STK_best_track_clusters->Fill(selected_track->getNhitX(), weight);

            std::vector<int> track_nHoles(2, 0);
            std::vector<double> track_slope(2, 0);
            std::vector<double> track_intercept(2, 0);
            std::vector<double> extr_BGO_top(2, 0);

            // Get the X and Y clusters
            ps_histos->h_STK_X_clusters_best_track->Fill(selected_track->getNhitX(), weight);
            ps_histos->h_STK_Y_clusters_best_track->Fill(selected_track->getNhitY(), weight);
            ps_histos->h_STK_clusters_best_track->Fill(selected_track->getNhitX(), selected_track->getNhitY(), weight);
            ps_histos->h_STK_X_clusters_vs_energy_best_track->Fill(evt_corr_energy_gev, selected_track->getNhitX(), weight);
            ps_histos->h_STK_Y_clusters_vs_energy_best_track->Fill(evt_corr_energy_gev, selected_track->getNhitY(), weight);
            
            track_points(
                selected_track,
                stkclusters,
                LadderToLayer,
                track_nHoles);

            // Get the X and Y holes
            ps_histos->h_STK_X_holes_best_track->Fill(track_nHoles[0]);
            ps_histos->h_STK_Y_holes_best_track->Fill(track_nHoles[1]);
            ps_histos->h_STK_X_holes_vs_energy_best_track->Fill(evt_corr_energy_gev, track_nHoles[0], weight);
            ps_histos->h_STK_Y_holes_vs_energy_best_track->Fill(evt_corr_energy_gev, track_nHoles[1], weight);
            
            // Find slope and intercept
            track_slope[0] = selected_track->getTrackParams().getSlopeX();
            track_slope[1] = selected_track->getTrackParams().getSlopeY();
            track_intercept[0] = selected_track->getTrackParams().getInterceptX();
            track_intercept[1] = selected_track->getTrackParams().getInterceptY();
            TVector3 trackDirection = (selected_track->getDirection()).Unit();

            // Extrapolate to the top of BGO
            for (int coord = 0; coord < 2; ++coord)
                extr_BGO_top[coord] = track_slope[coord] * BGO_TopZ + track_intercept[coord];

            // Evaluate distance between Top STK and BGO points
            double dxTop = extr_BGO_top[0] - bgoRecEntrance[0];
            double dyTop = extr_BGO_top[1] - bgoRecEntrance[1];
            double drTop = sqrt(pow(dxTop, 2) + pow(dyTop, 2));
            // Evaluate angular distance between STK track and BGO Rec track
            double dAngleTrackBgoRec = trackDirection.Angle(bgoRecDirection) * TMath::RadToDeg();
            
            ps_histos->h_STK_BGO_TOP_spatial_difference_best_track->Fill(drTop, weight);
            ps_histos->h_STK_BGO_TOP_spatial_X_difference_best_track->Fill(dxTop, weight);
            ps_histos->h_STK_BGO_TOP_spatial_Y_difference_best_track->Fill(dyTop, weight);
            ps_histos->h_STK_BGO_track_angular_difference_best_track->Fill(dAngleTrackBgoRec, weight);

            if (selected_track->getNhitX() == 3 && selected_track->getNhitY() == 3) {
                ps_histos->h_STK_BGO_TOP_spatial_difference_3_clusters_best_track->Fill(drTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_X_difference_3_clusters_best_track->Fill(dxTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_Y_difference_3_clusters_best_track->Fill(dyTop, weight);
                ps_histos->h_STK_BGO_track_angular_difference_3_clusters_best_track->Fill(dAngleTrackBgoRec, weight);
            }
            else if (selected_track->getNhitX() == 4 && selected_track->getNhitY() == 4) {
                ps_histos->h_STK_BGO_TOP_spatial_difference_4_clusters_best_track->Fill(drTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_X_difference_4_clusters_best_track->Fill(dxTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_Y_difference_4_clusters_best_track->Fill(dyTop, weight);
                ps_histos->h_STK_BGO_track_angular_difference_4_clusters_best_track->Fill(dAngleTrackBgoRec, weight);
            }
            else if (selected_track->getNhitX() == 5 && selected_track->getNhitY() == 5) {
                ps_histos->h_STK_BGO_TOP_spatial_difference_5_clusters_best_track->Fill(drTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_X_difference_5_clusters_best_track->Fill(dxTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_Y_difference_5_clusters_best_track->Fill(dyTop, weight);
                ps_histos->h_STK_BGO_track_angular_difference_5_clusters_best_track->Fill(dAngleTrackBgoRec, weight);
            }
            else if (selected_track->getNhitX() == 6 && selected_track->getNhitY() == 6) {
                ps_histos->h_STK_BGO_TOP_spatial_difference_6_clusters_best_track->Fill(drTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_X_difference_6_clusters_best_track->Fill(dxTop, weight);
                ps_histos->h_STK_BGO_TOP_spatial_Y_difference_6_clusters_best_track->Fill(dyTop, weight);
                ps_histos->h_STK_BGO_track_angular_difference_6_clusters_best_track->Fill(dAngleTrackBgoRec, weight);
            }

            ps_histos->h_STK_clusters_vs_angular_distance_X_best_track->Fill(selected_track->getNhitX(), dAngleTrackBgoRec, weight);
            ps_histos->h_STK_clusters_vs_angular_distance_Y_best_track->Fill(selected_track->getNhitY(), dAngleTrackBgoRec, weight);

            // Get the charge
            stk_charge_explore(stkclusters, selected_track, ps_histos);
        }
    }

std::tuple<bool, DmpStkTrack> get_best_track(
	const std::shared_ptr<DmpEvtBgoRec> bgorec,
    std::shared_ptr<DmpEvtHeader> evt_header, 
    const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtBgoHits> bgohits,
	const std::shared_ptr<TClonesArray> stkclusters,
	const std::shared_ptr<TClonesArray> stktracks,
    const double evt_energy,
    std::shared_ptr<config> _config) {
        
        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
        std::unique_ptr<DmpStkContainer> stkVault = std::make_unique<DmpStkContainer>();
        
        best_track event_best_track;

        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), (_config->GetCutsConfig()).bgo_layer_min_energy);
        stkVault->scanSTKHits(stkclusters);

        bool isTrackInteresting {false};

        if (check_trigger(evt_header)) {

            auto maxelayer_cut = maxElayer_cut(bgoVault->GetLayerEnergies(), (_config->GetCutsConfig()).bgo_max_energy_ratio, evt_energy);
            auto maxbarlayer_cut = maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
            auto bgotrack_cut = BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), (_config->GetCutsConfig()).bgo_shower_axis_delta);

            auto bgofiducial_cut = maxelayer_cut && maxbarlayer_cut && bgotrack_cut;

            if (bgofiducial_cut) {
                auto nbarlayer13_cut = nBarLayer13_cut(bgohits, bgoVault->GetSingleLayerBarNumber(13), evt_energy);
                if (nbarlayer13_cut) {
                    auto maxrms_cut = maxRms_cut(bgoVault->GetELayer(), bgoVault->GetRmsLayer(), evt_energy, (_config->GetCutsConfig()).bgo_shower_width);

                    if (maxrms_cut) {

                        auto trackselection_cut = track_selection_cut(
                            bgorec, 
                            bgoVault->GetBGOslope(), 
                            bgoVault->GetBGOintercept(), 
                            bgohits, 
                            stkclusters, 
                            stktracks, 
                            event_best_track,
                            (_config->GetCutsConfig()).STK_BGO_delta_position,
                            (_config->GetCutsConfig()).STK_BGO_delta_track,
                            (_config->GetEventDisplayConfig()).min_track_X_clusters,
                            (_config->GetEventDisplayConfig()).min_track_Y_clusters,
                            (_config->GetCutsConfig()).track_X_holes,
                            (_config->GetCutsConfig()).track_Y_holes);

                        if (trackselection_cut) {
                            if (event_best_track.myBestTrack.getNhitX() == (_config->GetEventDisplayConfig()).
                            track_X_clusters && event_best_track.myBestTrack.getNhitY() == (_config->GetEventDisplayConfig()).track_Y_clusters) {
                                isTrackInteresting = true;
                            }
                        }
                    }
                }
            }
        }
            
        return std::tuple<bool, DmpStkTrack> (isTrackInteresting, event_best_track.myBestTrack);
    }