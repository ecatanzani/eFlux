#ifndef CHARGE_H
#define CHARGE_H

#include "data_cuts.h"

#include "TH1D.h"
#include "TH2D.h"

extern void fillChargeHistos(
    TH1D &h_chargeX, 
    TH1D &h_chargeY,
    TH1D &h_charge,
    TH2D &h_charge2D,
    const best_track track,
    const std::shared_ptr<TClonesArray> stkclusters);

#endif