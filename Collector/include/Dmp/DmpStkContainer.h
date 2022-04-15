#ifndef DMPSTKCONTAINER_H
#define DMPSTKCONTAINER_H

#include <vector>
#include <memory>

#include "TVector3.h"
#include "TClonesArray.h"

#include "Dmp/DmpGeoStruct.h"

#include "DmpStkTrack.h"
#include "DmpStkSiCluster.h"

class DmpStkContainer
{
public:
	DmpStkContainer() : clusters_on_plane(DAMPE_stk_planes, 0),
						stkEcore1Rm (2, 0),
						nStkClu1Rm (2, 0)
						{
						}
	~DmpStkContainer(){};

	void scanSTKHits(
		const std::shared_ptr<TClonesArray> stkclusters,
		const std::shared_ptr<TClonesArray> stktracks,
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept);
	void Reset();
	const std::vector<int> GetNPlaneClusters();
	const std::vector<double> GetStkEcore1Rm();
	const std::vector<unsigned int> GetNStkClu1Rm();
	
private:
	void clusters_1_rm(
		const std::shared_ptr<TClonesArray> stkclusters,
		const std::shared_ptr<TClonesArray> stktracks,
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept);
	void link_ladders(std::vector<int> &LadderToLayer);
	void fill_BGO_vectors(
		TVector3 &bgoRecEntrance,
		TVector3 &bgoRecDirection,
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept);
	void get_track_points(
		const std::shared_ptr<TClonesArray> stkclusters,
		DmpStkTrack *track,
		const std::vector<int> LadderToLayer,
		std::vector<int> &track_nHoles);

	std::vector<int> clusters_on_plane;
	std::vector<double> stkEcore1Rm;
	std::vector<unsigned int> nStkClu1Rm;

};

#endif
