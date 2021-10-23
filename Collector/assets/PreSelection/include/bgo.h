#ifndef BGO_H
#define BGO_H

#include <memory>

#include "histos.h"

#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtHeader.h"

extern void bgofiducial_distributions(
    std::shared_ptr<DmpEvtBgoHits> bgohits,
    std::shared_ptr<DmpEvtBgoRec> bgorec,
    std::shared_ptr<DmpEvtHeader> evt_header,
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos);

#endif