#ifndef STK_H
#define STK_H

#include <memory>

#include "histos.h"
#include "config.h"

#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtHeader.h"

#include "TClonesArray.h"

#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtHeader.h"
#include "DmpEvtPsdHits.h"

#include "DmpStkTrack.h"

extern void stk_distributions(
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
    std::shared_ptr<config> cuts_config);

extern void BGO_vectors(
	TVector3 &bgoRecEntrance,
	TVector3 &bgoRecDirection,
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept);

extern void ladders(std::vector<int> &LadderToLayer);

extern void track_points(
	DmpStkTrack *track,
	const std::shared_ptr<TClonesArray> stkclusters,
	const std::vector<int> LadderToLayer,
	std::vector<int> &track_nHoles);

extern void track(
	const std::shared_ptr<DmpEvtBgoRec> bgorec,
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtBgoHits> bgohits,
	const std::shared_ptr<TClonesArray> stkclusters,
	const std::shared_ptr<TClonesArray> stktracks,
    const double evt_corr_energy_gev,
    std::shared_ptr<histos> ps_histos);

#endif