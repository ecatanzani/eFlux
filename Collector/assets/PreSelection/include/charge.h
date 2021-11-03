#ifndef CHARGE_H
#define CHARGE_H

#include <memory>

#include "histos.h"
#include "cuts.h"

#include "TClonesArray.h"

#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtHeader.h"
#include "DmpEvtPsdHits.h"

extern void charge_distributions(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header, 
    std::shared_ptr<TClonesArray> stkclusters, 
    std::shared_ptr<TClonesArray> stktracks,
    std::shared_ptr<DmpEvtPsdHits> psdhits,
    const double evt_energy, 
    const double evt_corr_energy,
    const double evt_energy_gev, 
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos);

extern void stk_charge(
    const std::shared_ptr<TClonesArray> stkclusters, 
    best_track &event_best_track, 
    std::shared_ptr<histos> ps_histos,
    const bool lastcut = false);

extern void stk_charge_explore(
    const std::shared_ptr<TClonesArray> stkclusters,
    DmpStkTrack* track,
    std::shared_ptr<histos> ps_histos);

extern void psd_charge(
	const std::vector<std::vector<double>> psdCluster_maxE,
    best_track &event_best_track,
    psd_cluster_match &clu_matching,
	std::shared_ptr<histos> ps_histos,
    const bool lastcut = false);

extern void psd_charge_explore(
	const std::vector<std::vector<double>> psdCluster_maxE,
    best_track &event_best_track,
    psd_cluster_match &clu_matching,
	std::shared_ptr<histos> ps_histos);

#endif