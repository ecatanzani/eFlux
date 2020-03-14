#include "acceptance_cuts.h"
#include "acceptance.h"

bool maxElater_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf &acceptance_cuts,
    const double bgoTotalE)
{
    // Get energy maximum along X and Y views
    double ELayer_max_XZ = 0;
    double ELayer_max_YZ = 0;

    for (int lIdx = 1; lIdx < DAMPE_bgo_nLayers; lIdx += 2)
    {
        auto lEgy = bgorec->GetELayer(lIdx);
        if (lEgy > ELayer_max_XZ)
            ELayer_max_XZ = lEgy;
    }

    for (int lIdx = 1; lIdx < DAMPE_bgo_nLayers; lIdx += 1)
    {
        auto lEgy = bgorec->GetELayer(lIdx);
        if (lEgy > ELayer_max_YZ)
            ELayer_max_YZ = lEgy;
    }

    bool passed_maxELayerTotalE_cut = true;
    double MaxELayer;
    if (ELayer_max_XZ > ELayer_max_YZ)
        MaxELayer = ELayer_max_XZ;
    else
        MaxELayer = ELayer_max_YZ;
    double rMaxELayerTotalE = MaxELayer / bgoTotalE;
    if (rMaxELayerTotalE > acceptance_cuts.energy_lRatio)
        passed_maxELayerTotalE_cut = false;

    return passed_maxELayerTotalE_cut;
}

bool maxBarLayer_cut(
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const int nBgoHits)
{
    bool passed_maxBarLayer_cut = true;
    std::vector<short> barNumberMaxEBarLay1_2_3(3, -1); // Bar number of maxE bar in layer 1, 2, 3
    std::vector<double> MaxEBarLay1_2_3(3, 0);          // E of maxE bar in layer 1, 2, 3

    for (int ihit = 0; ihit < nBgoHits; ++ihit)
    {
        auto hitE = (bgohits->fEnergy)[ihit];
        auto lay = bgohits->GetLayerID(ihit);
        if (lay == 1 || lay == 2 || lay == 3)
        {
            if (hitE > MaxEBarLay1_2_3[lay - 1])
            {
                auto iBar = ((bgohits->fGlobalBarID)[ihit] >> 6) & 0x1f;
                MaxEBarLay1_2_3[lay - 1] = hitE;
                barNumberMaxEBarLay1_2_3[lay - 1] = iBar;
            }
        }
    }

    for (int j = 0; j < 3; ++j)
        if (barNumberMaxEBarLay1_2_3[j] <= 0 || barNumberMaxEBarLay1_2_3[j] == 21)
            passed_maxBarLayer_cut = false;

    return passed_maxBarLayer_cut;
}

bool BGOTrackContainment_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf &acceptance_cuts,
    bool &passEvent)
{
    bool passed_bgo_containment_cut = false;
    double BGO_TopZ = 46;
    double BGO_BottomZ = 448;
    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);
    bgoRec_slope[1] = bgorec->GetSlopeXZ();
    bgoRec_slope[0] = bgorec->GetSlopeYZ();
    bgoRec_intercept[1] = bgorec->GetInterceptXZ();
    bgoRec_intercept[0] = bgorec->GetInterceptYZ();

    if ((bgoRec_slope[1] == 0 && bgoRec_intercept[1] == 0) ||
        (bgoRec_slope[0] == 0 && bgoRec_intercept[0] == 0))
        passEvent = false;

    double topX = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];
    double topY = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
    double bottomX = bgoRec_slope[1] * BGO_BottomZ + bgoRec_intercept[1];
    double bottomY = bgoRec_slope[0] * BGO_BottomZ + bgoRec_intercept[0];

    if (fabs(topX) < acceptance_cuts.shower_axis_delta && fabs(topY) < acceptance_cuts.shower_axis_delta && fabs(bottomX) < acceptance_cuts.shower_axis_delta && fabs(bottomY) < acceptance_cuts.shower_axis_delta)
        passed_bgo_containment_cut = true;

    return passed_bgo_containment_cut;
}

bool nBarLayer13_cut(
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const std::vector<short> &layerBarNumber,
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
    const std::vector<std::vector<short>> &layerBarNumber,
    const std::vector<double> &rmsLayer,
    const double bgoTotalE,
    const acceptance_conf &acceptance_cuts)
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

bool track_selection_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    TClonesArray* stkclusters,
    TClonesArray* stktracks,
    const acceptance_conf &acceptance_cuts)
{
    bool passed_track_selection_cut = false;
    std::vector<DmpStkTrack*> result = std::vector<DmpStkTrack*> ();
    TVector3 angle_BGO(bgorec->GetSlopeXZ(), bgorec->GetSlopeYZ(), 1);
    
    // Get BGO top impact point
    double BGO_TopZ = 46;
    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);
    bgoRec_slope[1] = bgorec->GetSlopeXZ();
    bgoRec_slope[0] = bgorec->GetSlopeYZ();
    bgoRec_intercept[1] = bgorec->GetInterceptXZ();
    bgoRec_intercept[0] = bgorec->GetInterceptYZ();
    double BGO_topX = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];
    double BGO_topY = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];

    // Loop on the tracks
    for (int trIdx = 0; trIdx < stktracks->GetLast(); ++trIdx)
    {
        // Track params...
        int clusterX_counter = 0;
        int clusterY_counter = 0;
        int nHoles_X = 0;
        int nHoles_Y = 0;
        TVector3 track_direction;
        TVector3 track_impact_point;
        double track_slopeX;
        double track_slopeY;
        double STK_topX;
        double STK_topY;

        // *********************

        auto track = static_cast<DmpStkTrack *>(stktracks->ConstructedAt(trIdx));

        for (int tPoint = 0; tPoint < track->GetNPoints(); ++tPoint)
        {
            track_direction = track->getDirection();
            track_impact_point = track->getImpactPoint();
            track_slopeX = track->getTrackParams().getSlopeX();
            track_slopeY = track->getTrackParams().getSlopeY();

            // --> Count the number of clusters on both X and Y
            auto clx = track->GetClusterX(tPoint, stkclusters);
            auto cly = track->GetClusterY(tPoint, stkclusters);
            if (clx)
                ++clusterX_counter;
            if (cly)
                ++clusterY_counter;

            // --> Count the number of holes on both X and Y
            if (tPoint > 0 && tPoint < (track->GetNPoints() - 1))
            {
                if (!(track->getHitMeasX(tPoint) > -9999.))
                    ++nHoles_X;
                if (!(track->getHitMeasY(tPoint) > -9999.))
                    ++nHoles_Y;
            }
        }

        STK_topX = track_slopeX * BGO_TopZ + track_impact_point.x();
        STK_topY = track_slopeY * BGO_TopZ + track_impact_point.y();

        // Filtering track

        // CUT 0 - MINUMUM NUMBER OF CLUSTER
        if (clusterX_counter >= acceptance_cuts.track_X_clusters && clusterY_counter >= acceptance_cuts.track_Y_clusters)
            if (nHoles_X <= acceptance_cuts.track_missingHit_X && nHoles_Y <= acceptance_cuts.track_missingHit_Y)
                if (track_direction.Angle(angle_BGO) <= acceptance_cuts.STK_BGO_delta_track)
                    if (fabs(STK_topX - BGO_topX) <= acceptance_cuts.STK_BGO_delta_position && fabs(STK_topY - BGO_topY) <= acceptance_cuts.STK_BGO_delta_position)
                        result.push_back(track);
    }

    if (result.size())
        passed_track_selection_cut = true;

    return passed_track_selection_cut;
}