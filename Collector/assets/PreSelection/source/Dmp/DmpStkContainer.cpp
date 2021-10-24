#include "Dmp/DmpStkContainer.h"

void DmpStkContainer::scanSTKHits(const std::shared_ptr<TClonesArray> stkclusters)
{
    for (int idx=0; idx<stkclusters->GetEntries(); ++idx)
        ++clusters_on_plane[static_cast<DmpStkSiCluster*>(stkclusters->ConstructedAt(idx))->getPlane()];
}

const std::vector<int> DmpStkContainer::GetNPlaneClusters()
{
    return clusters_on_plane;
}