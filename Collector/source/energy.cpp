#include "energy.h"

energy::energy(
	const double input_simuEnergy,
	const double input_bgoTotalE_raw,
	const double input_bgoTotalE_corr,
	const double input_simu_energy_w,
	const double input_corr_energy_w)
{
	simuEnergy = input_simuEnergy;
	bgoTotalE_raw = input_bgoTotalE_raw;
	bgoTotalE_corr = input_bgoTotalE_corr;
	simu_energy_w = input_simu_energy_w;
	corr_energy_w = input_corr_energy_w;
}

void energy::Set(
	const double input_simuEnergy,
	const double input_bgoTotalE_raw,
	const double input_bgoTotalE_corr,
	const double input_simu_energy_w,
	const double input_corr_energy_w)
{
	simuEnergy = input_simuEnergy;
	bgoTotalE_raw = input_bgoTotalE_raw;
	bgoTotalE_corr = input_bgoTotalE_corr;
	simu_energy_w = input_simu_energy_w;
	corr_energy_w = input_corr_energy_w;
}

void energy::SetEnergies(
	const double input_simuEnergy,
	const double input_bgoTotalE_raw,
	const double input_bgoTotalE_corr)
{
	simuEnergy = input_simuEnergy;
	bgoTotalE_raw = input_bgoTotalE_raw;
	bgoTotalE_corr = input_bgoTotalE_corr;
}

void energy::Reset()
{
	simuEnergy = -999;
	bgoTotalE_raw = -999;
	bgoTotalE_corr = -999;
	simu_energy_w = -999;
	corr_energy_w = -999;
}

void energy::SetSimuEnergy(const double input_simuEnergy)
{
	simuEnergy = input_simuEnergy;
}

void energy::SetRawEnergy(const double input_bgoTotalE_raw)
{
	bgoTotalE_raw = input_bgoTotalE_raw;
}

void energy::SetCorrEnergy(const double input_bgoTotalE_corr)
{
	bgoTotalE_corr = input_bgoTotalE_corr;
}

void energy::SetSimuEnergyWeight(const double input_energy_w)
{
	simu_energy_w = input_energy_w;
}

void energy::SetCorrEnergyWeight(const double input_energy_w)
{
	corr_energy_w = input_energy_w;
}

const double energy::GetSimuEnergy()
{
	return simuEnergy;
}

const double energy::GetRawEnergy()
{
	return bgoTotalE_raw;
}

const double energy::GetCorrEnergy()
{
	return bgoTotalE_corr;
}

const double energy::GetSimuWeight()
{
	return simu_energy_w;
}

const double energy::GetCorrWeight()
{
	return corr_energy_w;
}