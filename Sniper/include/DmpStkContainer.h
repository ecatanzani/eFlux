#ifndef DMPSTKCONTAINER_H
#define DMPSTKCONTAINER_H

#include <vector>
#include <memory>

#include "TClonesArray.h"

#include "DAMPE_geo_structure.h"

#include "DmpStkSiCluster.h"

class DmpStkContainer
{
public:
	DmpStkContainer() : clusters_on_plane(DAMPE_stk_planes, 0)
						{
						}
	~DmpStkContainer(){};

	void scanSTKHits(const std::shared_ptr<TClonesArray> stkclusters);
	const std::vector<int> GetNPlaneClusters();
	
private:
	std::vector<int> clusters_on_plane;

};

#endif
