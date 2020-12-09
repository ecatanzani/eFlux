#ifndef MC_H
#define MC_H

#include <string>
#include <vector>

#include "TFile.h"
#include "TVector3.h"

#include "energy.h"
#include "anyoption.h"
#include "DmpNudContainer.h"
#include "DmpBgoContainer.h"
#include "DmpSimuTrajectory.h"
#include "DmpFilterContainer.h"

struct _tmp_filter
{
	filter_output output;
	best_track evt_best_track;
	psd_charge evt_psd_charge;
	stk_charge evt_stk_charge;
	bgo_classifiers evt_bgo_classifier;
	trigger_info evt_trigger_info;
};

struct _tmp_bgo
{
	std::vector<double> layer_energies;
	std::vector<std::vector<double>> layer_bar_energies;
	std::vector<double> slope;
	std::vector<double> intercept;
	TVector3 trajectory2D;
	double sumrms;
	std::vector<double> sumrms_layer;
	std::vector<double> energy_fraction_layer;
	double energy_fraction_last_layer;
	double energy_fraction_13th_layer;
	int last_energy_layer;
	int hits;
	std::vector<double> energy_1mr;
	std::vector<double> energy_2mr;
	std::vector<double> energy_3mr;
	std::vector<double> energy_5mr;
};

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
	std::vector<double> thuthtrajectory_x;
	std::vector<double> thuthtrajectory_y;
	std::vector<double> thuthtrajectory_z;
	std::vector<double> truthtrajectory_energy;
	std::vector<double> thuthtrajectory_start_x;
	std::vector<double> thuthtrajectory_start_y;
	std::vector<double> thuthtrajectory_start_z;
	std::vector<double> thuthtrajectory_stop_x;
	std::vector<double> thuthtrajectory_stop_y;
	std::vector<double> thuthtrajectory_stop_z;
	std::vector<double> truthtrajectory_trackID;
	std::vector<double> truthtrajectory_parentID;
	std::vector<double> truthtrajectory_charge;
	std::vector<double> truthtrajectory_PDG;
	std::vector<double> truthtrajectory_stop_index;
};

struct _tmp_energy
{
	double simu;
	double simu_w;
	double raw;
	double correct;
	double correct_w;
};

struct _tmp_nud
{
	std::vector<int> adc;
	int total_adc;
	int max_adc;
	int max_channel_ID;
};

extern void mcCore(
	const std::string inputPath,
	const std::string outputPath,
	const bool _VERBOSE,
	const bool pedantic,
	AnyOption &opt,
	const std::string wd);

extern void mcLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool verbose,
	const std::string wd);

extern std::shared_ptr<_tmp_filter> fillFilterTmpStruct(DmpFilterContainer &filter);

extern std::shared_ptr<_tmp_simu> fillSimuTmpStruct(
	const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
	const std::shared_ptr<TClonesArray> simu_trajectories);

extern std::shared_ptr<_tmp_bgo> fillBGOTmpStruct(DmpBgoContainer &bgoVault);

extern std::shared_ptr<_tmp_energy> fillEnergyTmpStruct(energy &evt_energy);

extern std::shared_ptr<_tmp_nud> fillNUDTmpStruct(DmpNudContainer &nudVault);

#endif