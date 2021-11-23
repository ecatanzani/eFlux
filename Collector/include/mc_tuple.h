#ifndef MC_TUPLE_H
#define MC_TUPLE_H

#include "mc.h"
#include "tuple.h"
#include "tmpstruct.h"
#include "mc_tmpstruct.h"
#include "preselection.h"

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
		const std::vector<int> _stk_clusters_on_plane,
		const std::shared_ptr<_tmp_bgo> _bgo_res,
		const std::shared_ptr<_tmp_simu> _simu_res,
		const std::shared_ptr<_tmp_energy> _energy_res,
		const std::shared_ptr<_tmp_nud> _nud_res,
		const p_cuts &preselection_cuts);

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
	double energy_corr_w = -999;
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
	std::vector<double> simu_thruthtrajectory_x;
	std::vector<double> simu_thruthtrajectory_y;
	std::vector<double> simu_thruthtrajectory_z;
	std::vector<double> simu_thruthtrajectory_energy;
	std::vector<double> simu_thruthtrajectory_start_x;
	std::vector<double> simu_thruthtrajectory_start_y;
	std::vector<double> simu_thruthtrajectory_start_z;
	std::vector<double> simu_thruthtrajectory_stop_x;
	std::vector<double> simu_thruthtrajectory_stop_y;
	std::vector<double> simu_thruthtrajectory_stop_z;
	std::vector<double> simu_thruthtrajectory_trackID;
	std::vector<double> simu_thruthtrajectory_parentID;
	std::vector<double> simu_thruthtrajectory_charge;
	std::vector<double> simu_thruthtrajectory_PDG;
	std::vector<double> simu_thruthtrajectory_stop_index;

	// Filters
	bool evtfilter_geometric_before_trigger = false;
	bool evtfilter_trigger_check = false;
};

#endif