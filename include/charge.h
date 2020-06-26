#ifndef CHARGE_H
#define CHARGE_H

#include "data_cuts.h"

#include "TH1D.h"

extern void fillChargeHistos(
    TH1D &h_chargeX, 
    TH1D &h_chargeY, 
    const best_track track,
    const std::shared_ptr<TClonesArray> stkclusters);

#endif