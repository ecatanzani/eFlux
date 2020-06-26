#include "charge.h"

extern void fillChargeHistos(
    TH1D &h_chargeX, 
    TH1D &h_chargeY, 
    const best_track track,
    const std::shared_ptr<TClonesArray> stkclusters)
{
    double cluster_chargeX = -999;
    double cluster_chargeY = -999;
    auto track_correction = track.myBestTrack->getDirection().CosTheta();

    for(auto clIdx=0; clIdx<track.myBestTrack->GetNPoints(); ++clIdx)
    {
        auto cluster_x = track.myBestTrack->GetClusterX(clIdx, stkclusters.get());
        auto cluster_y = track.myBestTrack->GetClusterY(clIdx, stkclusters.get());
        if (cluster_x && !cluster_x->getPlane())
            h_chargeX.Fill(cluster_x->getEnergy()*track_correction);
        if (cluster_y && !cluster_y->getPlane())
            h_chargeY.Fill(cluster_y->getEnergy()*track_correction);
    }
}