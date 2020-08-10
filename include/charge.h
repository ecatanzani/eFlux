#ifndef CHARGE_H
#define CHARGE_H

#include "data_cuts.h"

#include "TH1D.h"
#include "TH2D.h"

template <typename InputChargeStruct> extern void fillChargeHistos(
    TH1D &h_chargeX,
	TH1D &h_chargeY,
	TH1D &h_charge,
	TH2D &h_charge2D,
    const InputChargeStruct charge_struct,
	const double energy_w = 1);

extern template void fillChargeHistos(
    TH1D &h_chargeX,
	TH1D &h_chargeY,
	TH1D &h_charge,
	TH2D &h_charge2D,
	const stk_charge extracted_stk_charge,
	const double energy_w = 1);

extern template void fillChargeHistos(
	TH1D &h_chargeX,
	TH1D &h_chargeY,
	TH1D &h_charge,
	TH2D &h_charge2D,
	const psd_charge extracted_psd_charge,
	const double energy_w = 1);

#endif