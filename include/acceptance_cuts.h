#ifndef ACCEPTANCE_CUTS_H
#define ACCEPTANCE_CUTS_H

#include <memory>
#include <vector>

#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"

#include "acceptance.h"
#include "data_cuts.h"

// DAMPESW includes
#include "DmpEvtSimuPrimaries.h"
#include "DmpEvtBgoRec.h"

extern bool geometric_cut(const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries);
extern bool geometric_top_cut(const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries);

// Analysis functions
extern void evaluateTopBottomPosition(
    const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    TH1D &h_BGOrec_topX_vs_realX,
    TH1D &h_BGOrec_topY_vs_realY,
    TH1D &h_real_slopeX,
    TH1D &h_real_slopeY,
    TH1D &h_BGOrec_slopeX,
    TH1D &h_BGOrec_slopeY,
    TH1D &h_real_interceptX,
    TH1D &h_real_interceptY,
    TH1D &h_BGOrec_interceptX,
    TH1D &h_BGOrec_interceptY,
    TH2D &h_real_topMap,
    TH2D &h_BGOreco_topMap,
    TH2D &h_real_bottomMap,
    TH2D &h_BGOreco_bottomMap);

extern bool checkBGOreco(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries);

extern void fillExternalMap(
    const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries, 
    TH2D &h_noBGOenergy_real_topMap);

#endif