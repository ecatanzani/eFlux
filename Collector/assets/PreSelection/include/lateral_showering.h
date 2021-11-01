#ifndef LATERAL_SHOWERING_H
#define LATERAL_SHOWERING_H

#include <memory>

#include "histos.h"

#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtHeader.h"
#include "DmpEvtSimuPrimaries.h"

#include "TClonesArray.h"

#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtHeader.h"
#include "DmpEvtPsdHits.h"

extern void lateral_showering_distributions(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header,
    const double evt_energy, 
    const double evt_corr_energy,
    const double evt_energy_gev, 
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos);

extern void lateral_showering_distributions_lastcut(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header, 
    std::shared_ptr<TClonesArray> stkclusters, 
    std::shared_ptr<TClonesArray> stktracks,
    std::shared_ptr<DmpEvtPsdHits> psdhits,
    std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
    const double evt_energy, 
    const double evt_corr_energy,
    const double evt_energy_gev, 
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos);

#endif