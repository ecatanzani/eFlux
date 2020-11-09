#ifndef MC_TUPLE_H
#define MC_TUPLE_H

#include "tuple.h"
#include "mc_histos.h"

#include <memory>
#include <string>

#include "TChain.h"
#include "TVector3.h"

class mc_tuple : public ntuple
{
public:
	mc_tuple(const std::shared_ptr<TChain> chain);
	~mc_tuple(){};
	void GetEntry(int idx);
	void InitHistos(const std::vector<float> binning);
	void WriteHistos(const std::string output_path);
private:
	void set_simu_address();
	void fill_histos();

	std::shared_ptr<mc_histos> _histos;
	// Energy
	double weight_shift;
	double simu_energy;
	double simu_energy_w;
	double corr_energy_w;
	TVector3 *simuPosition = nullptr;
	TVector3 *simuMomentum = nullptr;
	double simuSlopeX;
	double simuSlopeY;
	double simuInterceptX;
	double simuInterceptY;

	// Filters
	bool evtfilter_geometric_before_trigger;
	bool evtfilter_trigger_check;
};

#endif