#include "MC/mc_tmpstruct.h"

std::shared_ptr<_tmp_simu> fillSimuTmpStruct(
	const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
	const std::shared_ptr<TClonesArray> simu_trajectories)
{
	std::shared_ptr<_tmp_simu> _simu_res = std::make_shared<_tmp_simu>();
	_simu_res->position = TVector3(simu_primaries->pv_x, simu_primaries->pv_y, simu_primaries->pv_z);
	_simu_res->momentum = TVector3(simu_primaries->pvpart_px, simu_primaries->pvpart_py, simu_primaries->pvpart_pz);
	_simu_res->radius = simu_primaries->pv_r;
	_simu_res->theta = simu_primaries->pv_theta;
	_simu_res->phi = simu_primaries->pv_phi;
	_simu_res->flux_w = simu_primaries->pv_fluxw;
	_simu_res->n_particle = simu_primaries->npvpart;
	_simu_res->cos_x = simu_primaries->pvpart_cosx;
	_simu_res->cos_y = simu_primaries->pvpart_cosy;
	_simu_res->cos_z = simu_primaries->pvpart_cosz;
	_simu_res->charge = simu_primaries->pvpart_q;
	_simu_res->zenith = simu_primaries->pvpart_zenith;
	_simu_res->azimuth = simu_primaries->pvpart_azimuth;
	_simu_res->w = simu_primaries->pvpart_weight;
	_simu_res->PDG = simu_primaries->pvpart_pdg;
	_simu_res->geocut = simu_primaries->pvpart_geocut;

	auto n_simu_trajectories = simu_trajectories->GetEntries();

	_simu_res->thruthtrajectory_x.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_y.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_z.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_energy.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_start_x.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_start_y.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_start_z.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_stop_x.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_stop_y.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_stop_z.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_trackID.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_parentID.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_charge.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_PDG.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_stop_index.resize(n_simu_trajectories);

	for (int trIdx = 0; trIdx<n_simu_trajectories; ++trIdx)
	{
		auto simu_track = static_cast<DmpSimuTrajectory*>(simu_trajectories->ConstructedAt(trIdx));

		_simu_res->thruthtrajectory_x[trIdx] = simu_track->px;
		_simu_res->thruthtrajectory_y[trIdx] = simu_track->py;
		_simu_res->thruthtrajectory_z[trIdx] = simu_track->pz;
		_simu_res->thruthtrajectory_energy[trIdx] = simu_track->ekin;
		_simu_res->thruthtrajectory_start_x[trIdx] = simu_track->start_x;
		_simu_res->thruthtrajectory_start_y[trIdx] = simu_track->start_y;
		_simu_res->thruthtrajectory_start_z[trIdx] = simu_track->start_z;	
		_simu_res->thruthtrajectory_stop_x[trIdx] = simu_track->stop_x;
		_simu_res->thruthtrajectory_stop_y[trIdx] = simu_track->stop_y;
		_simu_res->thruthtrajectory_stop_z[trIdx] = simu_track->stop_z;
		_simu_res->thruthtrajectory_trackID[trIdx] = simu_track->trackID;
		_simu_res->thruthtrajectory_parentID[trIdx] = simu_track->parentID;
		_simu_res->thruthtrajectory_charge[trIdx] = simu_track->charge;
		_simu_res->thruthtrajectory_PDG[trIdx] = simu_track->pdg_id;
		_simu_res->thruthtrajectory_stop_index[trIdx] = simu_track->stop_index;
	}

	return _simu_res;
}

std::shared_ptr<_tmp_energy> fillEnergyTmpStruct(energy &evt_energy)
{
	std::shared_ptr<_tmp_energy> _energy_res = std::make_shared<_tmp_energy>();

	_energy_res->simu = evt_energy.GetSimuEnergy();
	_energy_res->simu_w = evt_energy.GetSimuWeight();
	_energy_res->raw = evt_energy.GetRawEnergy();
	_energy_res->correct = evt_energy.GetCorrEnergy();
	_energy_res->correct_w = evt_energy.GetCorrWeight();

	return _energy_res;
}