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
    TClonesArray *stkclusters,
    TClonesArray *stktracks,
    const acceptance_conf &acceptance_cuts)
{
    bool passed_track_selection_cut = false;
    std::vector<DmpStkTrack *> result = std::vector<DmpStkTrack *>();
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

bool xtrl_cut(
    const double sumRms,
    const std::vector<double> fracLayer,
    const acceptance_conf &acceptance_cuts)
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
    const acceptance_conf &acceptance_cuts)
{
    bool passed_psd_charge_cut = false;

    std::vector<short> vshort;
    std::vector<double> vdouble;
    std::vector<std::vector<short>> layerBarIndexPsd;     // Arrange PSD hit index
    std::vector<std::vector<short>> layerBarNumberPsd;    // Arrange PSD unique bar number
    std::vector<std::vector<double>> layerBarEnergyPsd;   // Arrange PSD hit energy
    std::vector<std::vector<short>> layerBarUsedPsd;      // mark PSD hits index used for clustering
    std::vector<std::vector<double>> psdClusterMaxECoord; // Arrange X and Y coordinates regarding the max energy release on PSD layer

    // Find max index and energy release for PSD hits
    std::vector<unsigned int> iMaxBarPsd(2, -1);
    std::vector<double> eMaxBarPsd(2, 0);
    bool first_hit = true;

    // Create vector matrix
    for (int nLayer = 0; nLayer < 2; ++nLayer)
    {
        layerBarIndexPsd.push_back(vshort);
        layerBarNumberPsd.push_back(vshort);
        layerBarEnergyPsd.push_back(vdouble);
        layerBarUsedPsd.push_back(vshort);
        psdClusterMaxECoord.push_back(vdouble);
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
            psdClusterMaxECoord[nLayer].push_back(actual_coordinate);
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

            /*
            psdClusterIndBeg[lay].push_back(indFirst);
            psdClusterLength[lay].push_back(cluSize);
            psdClusterIndMaxE[lay].push_back(imax);
            psdClusterMaxE[lay].push_back(eMaxBarPsd[lay]);
            psdClusterE[lay].push_back(cluE);
            psdClusterCoord[lay].push_back(thisCoord);
            psdClusterZ[lay].push_back(thisZ);
            */

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

    return passed_psd_charge_cut;
}