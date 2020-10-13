#include "charge.h"


template <typename InputChargeStruct> void fillChargeHistos(
    TH1D &h_chargeX,
	TH1D &h_chargeY,
	TH1D &h_charge,
	TH2D &h_charge2D,
    const InputChargeStruct charge_struct,
	const double energy_w)
{
	// Check charges
	if (charge_struct.chargeX != -999 && charge_struct.chargeY != -999)
	{
		// Compute mean charge
		auto mean_charge = 0.5 * (charge_struct.chargeX + charge_struct.chargeY);

		// Fill histos
		h_chargeX.Fill(charge_struct.chargeX, energy_w);
		h_chargeY.Fill(charge_struct.chargeY, energy_w);
		h_charge.Fill(mean_charge, energy_w);
		h_charge2D.Fill(charge_struct.chargeX, charge_struct.chargeY, energy_w);
	}
}

template void fillChargeHistos(
    TH1D &h_chargeX,
	TH1D &h_chargeY,
	TH1D &h_charge,
	TH2D &h_charge2D,
	const stk_charge extracted_stk_charge,
	const double energy_w = 1);

template void fillChargeHistos(
	TH1D &h_chargeX,
	TH1D &h_chargeY,
	TH1D &h_charge,
	TH2D &h_charge2D,
	const psd_charge extracted_psd_charge,
	const double energy_w = 1);