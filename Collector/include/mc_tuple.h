#ifndef MC_TUPLE_H
#define MC_TUPLE_H

#include "tuple.h"

#include "TVector3.h"

class mc_tuple : public ntuple
{
public:
	mc_tuple(const active_cuts &acuts)
	{
		init(acuts);
	};
	~mc_tuple(){};
	void Fill(
		const filter_output &output,
		const best_track &event_best_track,
		const double raw_energy,
		const double corr_energy,
		const double mc_corr_energy_w,
		const std::vector<double> &energy_release_layer,
		const std::vector<double> &bgoRec_slope,
		const std::vector<double> &bgoRec_intercept,
		const TVector3 &bgo_trajectory2D,
		const double sumRMS,
		const std::vector<double> &rms_layer,
		const std::vector<double> &bgo_fracLayer,
		const double lastFracLayer,
		const double frac_layer_13,
		const int last_bgo_layer,
		const int bgo_entries,
		const std::vector<double> &energy_1_moliere_radius,
		const std::vector<double> &energy_2_moliere_radius,
		const std::vector<double> &energy_3_moliere_radius,
		const std::vector<double> &energy_5_moliere_radius,
		const TVector3 &mc_position,
		const TVector3 &mc_momentum,
		const double mc_simu_energy,
		const double mc_simu_energy_w,
		const psd_charge &extracted_psd_charge,
		const stk_charge &extracted_stk_charge,
		const bgo_classifiers &classifier,
		const trigger_info &evt_trigger,
		const std::vector<double> &nud_adc,
		const double nud_total_adc,
		const double nud_max_adc,
		const int nud_max_channel_id);
	void Reset();

private:
	void init(const active_cuts &acuts);
	void branch_tree();
	void fill_filter_info(const filter_output &output);
	void fill_simu_info(
		const TVector3 mc_position,
		const TVector3 mc_momentum,
		const double mc_simu_energy,
		const double mc_corr_energy_w,
		const double mc_simu_energy_w);
	void reset_simu_info();

	// Energy
	double simu_energy = -999;
	double simu_energy_w = -999;
	double corr_energy_w = -999;
	TVector3 simuPosition{-999, -999, -999};
	TVector3 simuMomentum{-999, -999, -999};
	double simuSlopeX = -999;
	double simuSlopeY = -999;
	double simuInterceptX = -999;
	double simuInterceptY = -999;

	// Filters
	bool evtfilter_geometric_before_trigger = false;
	bool evtfilter_trigger_check = false;
};

#endif