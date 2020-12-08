#ifndef MC_TUPLE_H
#define MC_TUPLE_H

#include "mc.h"
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
		const std::shared_ptr<_tmp_filter> _filter_res,
		const std::shared_ptr<_tmp_bgo> _bgo_res,
		const std::shared_ptr<_tmp_simu> _simu_res,
		const std::shared_ptr<_tmp_energy> _energy_res,
		const std::shared_ptr<_tmp_nud> _nud_res);

	void Reset();

private:
	void init(const active_cuts &acuts);
	void branch_tree();
	void fill_filter_info(const filter_output &output);
	void fill_simu_info(
		const std::shared_ptr<_tmp_simu> _simu_res,
		const std::shared_ptr<_tmp_energy> _energy_res);
	void reset_simu_info();

	// Simu
	double simu_energy = -999;
	double simu_energy_w = -999;
	double corr_energy_w = -999;
	TVector3 simu_position{-999, -999, -999};
	TVector3 simu_momentum{-999, -999, -999};
	double simu_slope_x = -999;
	double simu_slope_y = -999;
	double simu_intercept_x = -999;
	double simu_intercept_y = -999;

	// Generation info
	double simu_radius = -999;
	double simu_theta = -999;
	double simu_phi = -999;
	double simu_flux_w = -999;
	int simu_n_particle = -999;
	double simu_cos_x = -999;
	double simu_cos_y = -999;
	double simu_cos_z = -999;
	double simu_charge = -999;
	double simu_zenith = -999;
	double simu_azimuth = -999;
	double simu_w = -999;
	long simu_PDG = -999;
	int simu_geocut = -999;

	// Simu Truth Trajectories
	double simu_thuthtrajectory_x = -999;
	double simu_thuthtrajectory_y = -999;
	double simu_thuthtrajectory_z = -999;
	double simu_truthtrajectory_energy = -999;
	double simu_thuthtrajectory_start_x = -999;
	double simu_thuthtrajectory_start_y = -999;
	double simu_thuthtrajectory_start_z = -999;
	double simu_thuthtrajectory_stop_x = -999;
	double simu_thuthtrajectory_stop_y = -999;
	double simu_thuthtrajectory_stop_z = -999;
	double simu_truthtrajectory_trackID = -999;
	double simu_truthtrajectory_parentID = -999;
	double simu_truthtrajectory_charge = -999;
	double simu_truthtrajectory_PDG = -999;
	double simu_truthtrajectory_stop_index = -999;

	// Filters
	bool evtfilter_geometric_before_trigger = false;
	bool evtfilter_trigger_check = false;
};

#endif