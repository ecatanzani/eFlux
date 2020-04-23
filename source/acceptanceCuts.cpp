#include "acceptance_cuts.h"
#include "acceptance.h"

bool checkBGOreco(const std::shared_ptr<DmpEvtBgoRec> bgorec)
{
    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);

    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

    if ((bgoRec_slope[0]==0 && bgoRec_intercept[0]==0) || (bgoRec_slope[1]==0 && bgoRec_intercept[1]==0))
        return false;
    else
        return true;
}

void fillExternalMap(
    const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries, 
    TH2D &h_noBGOenergy_real_topMap)
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

    h_noBGOenergy_real_topMap.Fill(actual_X, actual_Y);
}

bool geometric_cut(const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries)
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

    double actual_X = slope[0] * BGO_TopZ + intercept[0];
    double actual_Y = slope[1] * BGO_TopZ + intercept[1];

    if (fabs(actual_X) < BGO_SideXY && fabs(actual_Y) < BGO_SideXY)
        passed_geometric_cut = true;

#endif

    return passed_geometric_cut;
}

void evaluateTopPosition(
    const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    TH1D &h_BGOrec_topX_vs_realX,
    TH1D &h_BGOrec_topY_vs_realY,
    TH1D &h_real_slopeX,
    TH1D &h_real_slopeY,
    TH1D &h_BGOrec_slopeX,
    TH1D &h_BGOrec_slopeY,
    TH1D &h_real_interceptX,
    TH1D &h_real_interceptY,
    TH1D &h_BGOrec_interceptX,
    TH1D &h_BGOrec_interceptY,
    TH2D &h_real_topMap,
    TH2D &h_BGOreco_topMap)
{
    // Get the real simu position
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

    // Get the reco position
    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);

    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

    if ((bgoRec_slope[0]==0 && bgoRec_intercept[0]==0) || (bgoRec_slope[1]==0 && bgoRec_intercept[1]==0))
        return;

    double topZ = BGO_TopZ;
    double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
    double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

    /*
    double bottomZ = BGO_BottomZ;
    double bottomX = bgoRec_slope[0] * BGO_BottomZ + bgoRec_intercept[0];
    double bottomY = bgoRec_slope[1] * BGO_BottomZ + bgoRec_intercept[1];
    */

    // Fill slopes
    h_real_slopeX.Fill(slope[0]);
    h_real_slopeY.Fill(slope[1]);
    h_BGOrec_slopeX.Fill(bgoRec_slope[0]);
    h_BGOrec_slopeY.Fill(bgoRec_slope[1]);

    // Fill intercepts
    h_real_interceptX.Fill(intercept[0]);
    h_real_interceptY.Fill(intercept[1]);
    h_BGOrec_interceptX.Fill(bgoRec_intercept[0]);
    h_BGOrec_interceptY.Fill(bgoRec_intercept[1]);

    auto spreadX = actual_X - topX;
    auto spreadY = actual_Y - topY;

    // Fill spreads
    h_BGOrec_topX_vs_realX.Fill(spreadX);
    h_BGOrec_topY_vs_realY.Fill(spreadY);

    // Fill maps
    h_real_topMap.Fill(actual_X, actual_Y);
    h_BGOreco_topMap.Fill(topX, topY);
}

bool maxElayer_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf acceptance_cuts,
    const double bgoTotalE)
{
    bool passed_maxELayerTotalE_cut = true;

    int iMaxELayer = -1;  // Index of the layer corresponding to the max energy
    double MaxELayer = 0; // Value of the max energy

    // Found the max energy value and layer
    for (int idxLy = 0; idxLy < DAMPE_bgo_nLayers; ++idxLy)
    {
        auto layer_energy = static_cast<double>((bgorec->GetLayerEnergy())[idxLy]);
        if (layer_energy > MaxELayer)
        {
            MaxELayer = layer_energy;
            iMaxELayer = idxLy;
        }
    }

    auto rMaxELayerTotalE = MaxELayer / bgoTotalE;
    if (rMaxELayerTotalE > acceptance_cuts.energy_lRatio)
        passed_maxELayerTotalE_cut = false;

    return passed_maxELayerTotalE_cut;
}

void evaluateEnergyRatio(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf acceptance_cuts,
    const double bgoTotalE,
    TH1D &h_layer_energy_ratio)
{
    int iMaxELayer = -1;  // Index of the layer corresponding to the max energy
    double MaxELayer = 0; // Value of the max energy

    // Found the max energy value and layer
    for (int idxLy = 0; idxLy < DAMPE_bgo_nLayers; ++idxLy)
    {
        auto layer_energy = static_cast<double>((bgorec->GetLayerEnergy())[idxLy]);
        if (layer_energy > MaxELayer)
        {
            MaxELayer = layer_energy;
            iMaxELayer = idxLy;
        }
    }

    auto rMaxELayerTotalE = MaxELayer / bgoTotalE;
    h_layer_energy_ratio.Fill(rMaxELayerTotalE);
}

bool maxBarLayer_cut(
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

bool BGOTrackContainment_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf acceptance_cuts)
{
    bool passed_bgo_containment_cut = false;

    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);

    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

    if ((bgoRec_slope[0]==0 && bgoRec_intercept[0]==0) || (bgoRec_slope[1]==0 && bgoRec_intercept[1]==0))
        return passed_bgo_containment_cut;

    double topZ = BGO_TopZ;
    double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
    double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

    double bottomZ = BGO_BottomZ;
    double bottomX = bgoRec_slope[0] * BGO_BottomZ + bgoRec_intercept[0];
    double bottomY = bgoRec_slope[1] * BGO_BottomZ + bgoRec_intercept[1];

    if (
        fabs(topX) < acceptance_cuts.shower_axis_delta &&
        fabs(topY) < acceptance_cuts.shower_axis_delta &&
        fabs(bottomX) < acceptance_cuts.shower_axis_delta &&
        fabs(bottomY) < acceptance_cuts.shower_axis_delta)
        passed_bgo_containment_cut = true;
    
    return passed_bgo_containment_cut;
}

bool nBarLayer13_cut(
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

bool maxRms_cut(
    const std::vector<std::vector<short>> layerBarNumber,
    const std::vector<double> rmsLayer,
    const double bgoTotalE,
    const acceptance_conf acceptance_cuts)
{
    bool passed_maxRms_cut = false;
    auto max_rms = rmsLayer[0];
    auto eCut = bgoTotalE / 100.;

    for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
    {
        double layerTotEnergy = 0;
        std::accumulate(layerBarNumber[lIdx].begin(), layerBarNumber[lIdx].end(), layerTotEnergy);
        if (layerTotEnergy > eCut)
            if (rmsLayer[lIdx] > max_rms)
                max_rms = rmsLayer[lIdx];
    }
    if (max_rms < acceptance_cuts.max_rms_shower_width)
        passed_maxRms_cut = true;

    return passed_maxRms_cut;
}

inline void link_ladders(std::vector<int> &LadderToLayer)
{
    for (int ilad = 0; ilad < nSTKladders; ++ilad)
    {
        int iTRB = ilad / 24;
        int iladTRB = ilad % 24;
        int iPlane = 5 - iladTRB / 4; // 0-5 for each double layer
        // int isX = (iTRB / 2) % 2;     // measuring X, TRB2, 3, 6, 7
        int isY = (iTRB / 2 + 1) % 2; // measuring Y, TRB0, 1, 4, 5
        int iLay = iPlane * 2 + isY;  // 0-11 for each layer, starting from X (measuring X)
        LadderToLayer[ilad] = iLay;
    }
}

inline void fill_BGO_vectors(
    TVector3 &bgoRecEntrance,
    TVector3 &bgoRecDirection,
    const std::shared_ptr<DmpEvtBgoRec> bgorec)
{
    std::vector<double> bgoRec_slope(2, 0);
    std::vector<double> bgoRec_intercept(2, 0);

    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

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
    std::vector<unsigned int> &track_nHoles,
    best_track &event_best_track,
    const bool best_track = false)
{
    std::vector<int> prevHole(2 - 2);
    std::vector<int> firstLayer(2, -1);
    std::vector<int> lastLayer(2, -1);
    std::vector<int> lastPoint(2, -1);
    std::vector<unsigned int> track_nHoles_cont(2, 0);

    // Loop on track points to find last layer values
    for (int ip = track->GetNPoints() - 1; ip >= 0; --ip)
    {
        if (lastLayer[0] == -1)
        {
            if (track->getHitMeasX(ip) > -9999.)
            {
                lastPoint[0] = ip;
                DmpStkSiCluster *cluster = track->GetClusterX(ip, stkclusters.get());
                auto hardID = cluster->getLadderHardware();
                lastLayer[0] = LadderToLayer[hardID];
            }
        }
        if (lastLayer[1] == -1)
        {
            if (track->getHitMeasY(ip) > -9999.)
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
        if (track->getHitMeasX(ip) > -9999.)
        {
            DmpStkSiCluster *cluster = track->GetClusterX(ip, stkclusters.get());
            auto hardID = cluster->getLadderHardware();
            if (firstLayer[0] == -1)
                firstLayer[0] = LadderToLayer[hardID];
            continue;
        }
        if (firstLayer[0] != -1)
            ++track_nHoles[0];
        if (ip == prevHole[0] + 1)
            ++track_nHoles_cont[0];
        prevHole[0] = ip;
    }

    for (int ip = 0; ip <= lastPoint[1]; ++ip)
    {
        if (track->getHitMeasY(ip) > -9999.)
        {
            DmpStkSiCluster *cluster = track->GetClusterY(ip, stkclusters.get());
            auto hardID = cluster->getLadderHardware();
            if (firstLayer[1] == -1)
                firstLayer[1] = LadderToLayer[hardID];
            continue;
        }
        if (firstLayer[1] != -1)
            ++track_nHoles[1];
        if (ip == prevHole[1] + 1)
            ++track_nHoles_cont[1];
        prevHole[1] = ip;
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

bool track_selection_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const std::shared_ptr<TClonesArray> stkclusters,
    const std::shared_ptr<TClonesArray> stktracks,
    const acceptance_conf acceptance_cuts,
    best_track &event_best_track)
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
        bgorec);

    // Loop on the tracks
    for (int trIdx = 0; trIdx < stktracks->GetLast() + 1; ++trIdx)
    {
        std::vector<unsigned int> track_nHoles(2, 0);
        std::vector<double> track_slope(2, 0);
        std::vector<double> track_intercept(2, 0);
        std::vector<double> extr_BGO_top(2, 0);

        // *********************

        // Get the track
        auto track = static_cast<DmpStkTrack *>(stktracks->ConstructedAt(trIdx));

        // Reject tracks with not enough X and Y clusters
        if (track->getNhitX() < acceptance_cuts.track_X_clusters || track->getNhitY() < acceptance_cuts.track_Y_clusters)
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

        if (drTop > acceptance_cuts.STK_BGO_delta_position)
            continue;
        if (dAngleTrackBgoRec > acceptance_cuts.STK_BGO_delta_track)
            continue;

        selectedTracks.push_back(track);
    }

    // Sort selected tracks vector
    DmpStkTrackHelper tHelper(stktracks.get(), true, bgorec.get(), bgohits.get());
    tHelper.MergeSort(selectedTracks, &DmpStkTrackHelper::TracksCompare);

    if (selectedTracks.size() > 0)
    {
        DmpStkTrack *selected_track = static_cast<DmpStkTrack *>(selectedTracks[0]);
        std::vector<unsigned int> track_nHoles(2, 0);

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

        passed_track_selection_cut = true;
    }

    return passed_track_selection_cut;
}

bool xtrl_cut(
    const double sumRms,
    const std::vector<double> fracLayer,
    const acceptance_conf acceptance_cuts)
{
    bool passed_xtrl_cut = false;
    unsigned int last_layer_idx = 0;

    // Find the last layer with an energy release
    for (auto it = fracLayer.begin(); it != fracLayer.end(); ++it)
        if (*it)
        {
            auto index = std::distance(fracLayer.begin(), it);
            last_layer_idx = index;
        }

    // Build XTRL
    double xtrl = 0.1251e-8 * pow(sumRms, 4) * fracLayer[last_layer_idx];

    // Filter XTRL
    if (xtrl < acceptance_cuts.xtrl)
        passed_xtrl_cut = true;

    return passed_xtrl_cut;
}

bool psd_charge_cut(
    const std::shared_ptr<DmpEvtPsdHits> psdhits,
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf acceptance_cuts,
    const best_track event_best_track)
{
    bool passed_psd_charge_cut = false;

    // Getting BGO Reco slope and intercept (used for cluster matching)
    std::vector<double> bgoRec_slope(2, -999);
    std::vector<double> bgoRec_intercept(2, -999);
    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

    std::vector<short> vshort;
    std::vector<double> vdouble;
    std::vector<std::vector<short>> layerBarIndexPsd;   // Arrange PSD hit index
    std::vector<std::vector<short>> layerBarNumberPsd;  // Arrange PSD unique bar number
    std::vector<std::vector<double>> layerBarEnergyPsd; // Arrange PSD hit energy
    std::vector<std::vector<short>> layerBarUsedPsd;    // mark PSD hits index used for clustering

    // PSD clusters
    std::vector<std::vector<short>> psdCluster_idxBeg;
    std::vector<std::vector<short>> psdCluster_length;
    std::vector<std::vector<short>> psdCluster_idxMaxE;
    std::vector<std::vector<double>> psdCluster_E;
    std::vector<std::vector<double>> psdCluster_maxE;
    std::vector<std::vector<double>> psdCluster_maxEcoordinate; // Arrange X and Y coordinates regarding the max energy release on PSD layer
    std::vector<std::vector<double>> psdCluster_coordinate;     // weighted if more than 1 strip
    std::vector<std::vector<double>> psdCluster_Z;              // weighted if more than 1 strip

    // Find max index and energy release for PSD hits
    std::vector<unsigned int> iMaxBarPsd(2, -1);
    std::vector<double> eMaxBarPsd(2, 0);
    bool first_hit = true;

    // Create vector matrix
    for (int nLayer = 0; nLayer < 2; ++nLayer)
    {
        // PSD clustering...
        layerBarIndexPsd.push_back(vshort);
        layerBarNumberPsd.push_back(vshort);
        layerBarEnergyPsd.push_back(vdouble);
        layerBarUsedPsd.push_back(vshort);

        // PSD custers...
        psdCluster_idxBeg.push_back(vshort);
        psdCluster_length.push_back(vshort);
        psdCluster_idxMaxE.push_back(vshort);
        psdCluster_E.push_back(vdouble);
        psdCluster_maxE.push_back(vdouble);
        psdCluster_maxEcoordinate.push_back(vdouble);
        psdCluster_coordinate.push_back(vdouble);
        psdCluster_Z.push_back(vdouble);
    }

    // Loop on PSD hits
    int nPSD_tot_entries = psdhits->GetHittedBarNumber();
    for (int ihit = 0; ihit < nPSD_tot_entries; ++ihit)
    {
        double hitE = (psdhits->fEnergy)[ihit];
        if (hitE > 0.5)
        {
            short layerID = psdhits->GetLayerID(ihit);
            int iBar = ((psdhits->fGlobalBarID)[ihit] >> 6);
            layerBarUsedPsd[layerID].push_back(0);
            layerBarIndexPsd[layerID].push_back(ihit);
            layerBarNumberPsd[layerID].push_back(iBar);
            layerBarEnergyPsd[layerID].push_back(hitE);
            if (first_hit)
            {
                eMaxBarPsd[layerID] = hitE;
                iMaxBarPsd[layerID] = layerBarIndexPsd[layerID].size() - 1;
                first_hit = false;
            }
            else
            {
                if (hitE > eMaxBarPsd[layerID])
                {
                    eMaxBarPsd[layerID] = hitE;
                    iMaxBarPsd[layerID] = layerBarIndexPsd[layerID].size() - 1;
                }
            }
        }
    }

    /*
        cluster finding : starting from maxE, seed + the highest E neigboring bar
        coordinate: weighted average of two highest neigboring bars
    */

    for (int nLayer = 0; nLayer < 2; ++nLayer)
    {
        while (eMaxBarPsd[nLayer] > acceptance_cuts.PSD_bar_min_energy_release)
        {
            int psdhits_maxIdx = layerBarIndexPsd[nLayer][iMaxBarPsd[nLayer]];
            int psdbar_maxIdx = layerBarNumberPsd[nLayer][iMaxBarPsd[nLayer]];
            double cluster_energy = eMaxBarPsd[nLayer];
            double actual_coordinate = nLayer % 2 ? psdhits->GetHitX(psdhits_maxIdx) : psdhits->GetHitY(psdhits_maxIdx);
            psdCluster_maxEcoordinate[nLayer].push_back(actual_coordinate);
            double actualZ = psdhits->GetHitZ(psdhits_maxIdx);
            layerBarUsedPsd[nLayer][iMaxBarPsd[nLayer]] = 1;

            int cluster_closest_barIdx = -1;
            double energy_leftMaxBar = 0;
            double energy_rightMaxBar = 0;
            int cluster_fIdx = psdhits_maxIdx; // First index of the cluster
            int cluster_size = 1;              // Size of the cluster

            // Complete the cluster building procedure with the closest bars to the one with the maximum energy release
            if (iMaxBarPsd[nLayer] - 1 > 0) // Check the previous bar
                if ((layerBarUsedPsd[nLayer][iMaxBarPsd[nLayer] - 1] == 0) && (psdbar_maxIdx == layerBarNumberPsd[nLayer][iMaxBarPsd[nLayer] - 1] + 1))
                    energy_leftMaxBar = layerBarEnergyPsd[nLayer][iMaxBarPsd[nLayer] - 1];
            if ((iMaxBarPsd[nLayer] + 1) < layerBarEnergyPsd[nLayer].size()) // Check the next bar
                if ((layerBarUsedPsd[nLayer][iMaxBarPsd[nLayer] + 1] == 0) && (psdbar_maxIdx == layerBarNumberPsd[nLayer][iMaxBarPsd[nLayer] + 1] - 1))
                    energy_rightMaxBar = layerBarEnergyPsd[nLayer][iMaxBarPsd[nLayer] + 1];

            /* 
                FINAL CLUSTER BUILING:
                Now we have the central bar with the highest energy release, the left and the right one. 
                One cluster is composed by 2 bars. We have to choose, between the 2 lateral bars, that one with the highest energy release
            */
            if (energy_leftMaxBar > energy_rightMaxBar)
            {
                cluster_closest_barIdx = psdhits_maxIdx - 1;
                cluster_energy += energy_leftMaxBar;
                ++cluster_size;
                cluster_fIdx = psdhits_maxIdx - 1;
                layerBarUsedPsd[nLayer][iMaxBarPsd[nLayer] - 1] = 1;
            }
            //else if (energy_leftMaxBar < energy_rightMaxBar)
            else
            {
                cluster_closest_barIdx = psdhits_maxIdx + 1;
                cluster_energy += energy_rightMaxBar;
                ++cluster_size;
                layerBarUsedPsd[nLayer][iMaxBarPsd[nLayer] + 1] = 1;
            }

            if (cluster_closest_barIdx > -1)
            {
                double cluster_closest_bar_coordinate = iMaxBarPsd[nLayer] % 2 ? psdhits->GetHitX(cluster_closest_barIdx) : psdhits->GetHitY(cluster_closest_barIdx);
                double cluster_closest_energy = (psdhits->fEnergy)[cluster_closest_barIdx];
                double cluster_closestZ = psdhits->GetHitZ(cluster_closest_barIdx);
                actual_coordinate = (actual_coordinate * eMaxBarPsd[nLayer] + cluster_closest_bar_coordinate * cluster_closest_energy) / (eMaxBarPsd[nLayer] + cluster_closest_energy);
                actualZ = (actualZ * eMaxBarPsd[nLayer] + cluster_closestZ * cluster_closest_energy) / (eMaxBarPsd[nLayer] + cluster_closest_energy);
            }

            psdCluster_idxBeg[nLayer].push_back(cluster_fIdx);
            psdCluster_length[nLayer].push_back(cluster_size);
            psdCluster_idxMaxE[nLayer].push_back(psdhits_maxIdx);
            psdCluster_maxE[nLayer].push_back(eMaxBarPsd[nLayer]);
            psdCluster_E[nLayer].push_back(cluster_energy);
            psdCluster_coordinate[nLayer].push_back(actual_coordinate);
            psdCluster_Z[nLayer].push_back(actualZ);

            // Search for the next maximum energy release bar
            eMaxBarPsd[nLayer] = 0; // Reset maximum energy release
            for (unsigned int barIdx = 0; barIdx < layerBarEnergyPsd[nLayer].size(); ++barIdx)
            {
                // Check for the next unused PSD bar with a non negligible energy release
                if (layerBarUsedPsd[nLayer][barIdx] == 0 && layerBarEnergyPsd[nLayer][barIdx] > eMaxBarPsd[nLayer])
                {
                    iMaxBarPsd[nLayer] = barIdx;
                    eMaxBarPsd[nLayer] = layerBarEnergyPsd[nLayer][barIdx];
                }
            }
        }
    }

    // Get the number of clusters on both X and Y
    auto nPsdClusters = psdCluster_idxBeg[0].size() + psdCluster_idxBeg[1].size();
    auto nPsdClustersX = psdCluster_idxBeg[0].size();
    auto nPsdClustersY = psdCluster_idxBeg[1].size();

    // PSD clusters matching

    std::vector<int> icloPsdClu(2, -1); //0 is Y, 0 is X
    std::vector<double> dxCloPsdClu(2, 9999);
    std::vector<int> icloPsdCluMaxHit(2, -1);
    std::vector<double> dxCloPsdCluMaxHit(2, 9999);
    std::vector<int> icloPsdClu_bgoRec(2, -1);
    std::vector<double> dxCloPsdClu_bgoRec(2, 9999);
    std::vector<int> icloPsdClu_track(2, -1);
    std::vector<double> dxCloPsdClu_track(2, 9999);
    std::vector<int> icloPsdClu2_track(2, -1);
    std::vector<double> dxCloPsdClu2_track(2, 9999);

    for (int nLayer = 0; nLayer < 2; ++nLayer)
    {
        for (unsigned int iclu = 0; iclu < psdCluster_idxBeg[nLayer].size(); ++iclu)
        {
            bool IsMeasuringX = nLayer % 2;
            double hitZ = psdCluster_Z[nLayer][iclu];
            double thisCoord = psdCluster_maxEcoordinate[nLayer][iclu];

            // Get distance between actual coordinate and BGO rec coordinate
            double projCoord_bgoRec = IsMeasuringX ? bgoRec_slope[0] * hitZ + bgoRec_intercept[0] : bgoRec_slope[1] * hitZ + bgoRec_intercept[1];
            double dX_bgoRec = thisCoord - projCoord_bgoRec;
            if (fabs(dX_bgoRec) < fabs(dxCloPsdClu_bgoRec[nLayer]))
            {
                dxCloPsdClu_bgoRec[nLayer] = dX_bgoRec;
                icloPsdClu_bgoRec[nLayer] = iclu;
            }

            // Get distance between actual coordinate and best track coordinate
            double projCoord_track = IsMeasuringX ? event_best_track.track_slope[0] * hitZ + event_best_track.track_intercept[0] : event_best_track.track_slope[1] * hitZ + event_best_track.track_intercept[1];
            double dX_track = thisCoord - projCoord_track;

            if (fabs(dX_track) < fabs(dxCloPsdClu_track[nLayer]))
            {
                dxCloPsdClu_track[nLayer] = dX_track;
                icloPsdClu_track[nLayer] = iclu;
            }
        }
    }

    passed_psd_charge_cut = (fabs(dxCloPsdClu_track[0]) < acceptance_cuts.STK_PSD_delta_position && fabs(dxCloPsdClu_track[1]) < acceptance_cuts.STK_PSD_delta_position) ? true : false;
    return passed_psd_charge_cut;
}
