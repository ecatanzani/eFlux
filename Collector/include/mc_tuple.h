#ifndef MC_TUPLE_H
#define MC_TUPLE_H

#include "tuple.h"

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
		const std::vector<double> &energy_release_layer,
		const std::vector<double> &bgoRec_slope,
		const std::vector<double> &bgoRec_intercept,
		const double sumRMS,
		const std::vector<double> rms_layer,
		const std::vector<double> bgo_fracLayer,
		const double lastFracLayer,
		const double frac_layer_13,
		const int last_bgo_layer,
		const int bgo_entries,
		const std::vector<double> energy_1_moliere_radius,
		const std::vector<double> energy_2_moliere_radius,
		const std::vector<double> energy_3_moliere_radius,
		const std::vector<double> energy_5_moliere_radius,
		const psd_charge &extracted_psd_charge,
		const stk_charge &extracted_stk_charge,
		const bgo_classifiers &classifier,
		const trigger_info &evt_trigger,
		const std::vector<double> adc,
		const double total_adc,
		const double max_adc,
		const int max_channel_id);
	void Reset();
	
private:
	void init(const active_cuts &acuts);
	void branch_tree();
	void fill_filter_info(const filter_output &output);

	// Filters
	bool evtfilter_geometric_before_trigger = false;
	bool evtfilter_trigger_check = false;
};

#endif