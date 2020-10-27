#ifndef ENERGY_H
#define ENERGY_H

class energy
{
public:
	energy(){};
	energy(
		const double input_simuEnergy,
		const double input_bgoTotalE_raw,
		const double input_bgoTotalE_corr,
		const double input_energy_w);
	~energy(){};
	void Set(
		const double input_simuEnergy,
		const double input_bgoTotalE_raw,
		const double input_bgoTotalE_corr,
		const double input_energy_w);
	void SetEnergies(
		const double input_simuEnergy,
		const double input_bgoTotalE_raw,
		const double input_bgoTotalE_corr);
	void Reset();
	void SetSimuEnergy(const double input_simuEnergy);
	void SetRawEnergy(const double input_bgoTotalE_raw);
	void SetCorrEnergy(const double input_bgoTotalE_corr);
	void SetEnergyWeight(const double input_energy_w);
	const double GetSimuEnergy();
	const double GetRawEnergy();
	const double GetCorrEnergy();
	const double GetWeight();

private:
	double simuEnergy = -999;
	double bgoTotalE_raw = -999;
	double bgoTotalE_corr = -999;
	double energy_w = -999;
};

#endif