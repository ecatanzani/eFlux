#include "charge.h"

void fillChargeHistos(
    TH1D &h_chargeX, 
    TH1D &h_chargeY,
    TH1D &h_charge,
    TH2D &h_charge2D,
    const best_track track,
    const std::shared_ptr<TClonesArray> stkclusters)
{
    double cluster_chargeX = -999;
    double cluster_chargeY = -999;
    auto btrack = track.myBestTrack;
    
    // Charge correction 
    auto track_correction = btrack.getDirection().CosTheta();

    // Compute charges
    for(auto clIdx=0; clIdx<track.n_points; ++clIdx)
    {
        auto cluster_x = btrack.GetClusterX(clIdx, stkclusters.get());
        auto cluster_y = btrack.GetClusterY(clIdx, stkclusters.get());
        if (cluster_x && !cluster_x->getPlane())
            cluster_chargeX = cluster_x->getEnergy()*track_correction;
        if (cluster_y && !cluster_y->getPlane())
            cluster_chargeY = cluster_y->getEnergy()*track_correction;
    }

    // Check charges
    if (cluster_chargeX != -999 && cluster_chargeY != -999)
    {
        // Compute mean charge
        auto mean_charge = 0.5*(cluster_chargeX + cluster_chargeY);
        
        // Fill histos
        h_chargeX.Fill(cluster_chargeX);
        h_chargeY.Fill(cluster_chargeY);
        h_charge.Fill(mean_charge);
        h_charge2D.Fill(cluster_chargeX, cluster_chargeY);
    }

}