#ifndef MC_TMPSTRUCT_H
#define MC_TMPSTRUCT_H

#include <vector>
#include <memory>

#include "TVector3.h"
#include "TClonesArray.h"

#include "DmpSimuTrajectory.h"
#include "DmpEvtSimuPrimaries.h"

#include "energy.h"

struct _tmp_simu
{
	TVector3 position;
	TVector3 momentum;
	double radius;
	double theta;
	double phi;
	double flux_w;
	int n_particle;
	double cos_x;
	double cos_y;
	double cos_z;
	double charge;
	double zenith;
	double azimuth;
	double w;
	long PDG;
	int geocut;
	std::vector<double> thruthtrajectory_x;
	std::vector<double> thruthtrajectory_y;
	std::vector<double> thruthtrajectory_z;
	std::vector<double> thruthtrajectory_energy;
	std::vector<double> thruthtrajectory_start_x;
	std::vector<double> thruthtrajectory_start_y;
	std::vector<double> thruthtrajectory_start_z;
	std::vector<double> thruthtrajectory_stop_x;
	std::vector<double> thruthtrajectory_stop_y;
	std::vector<double> thruthtrajectory_stop_z;
	std::vector<double> thruthtrajectory_trackID;
	std::vector<double> thruthtrajectory_parentID;
	std::vector<double> thruthtrajectory_charge;
	std::vector<double> thruthtrajectory_PDG;
	std::vector<double> thruthtrajectory_stop_index;
};

struct _tmp_energy
{
	double simu;
	double simu_w;
	double raw;
	double correct;
	double correct_w;
};

extern std::shared_ptr<_tmp_simu> fillSimuTmpStruct(
	const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
	const std::shared_ptr<TClonesArray> simu_trajectories);
extern std::shared_ptr<_tmp_energy> fillEnergyTmpStruct(energy &evt_energy);

#endif