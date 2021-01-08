#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

struct cuts_conf
{
	double min_event_energy = -999;
	double max_event_energy = -999;
	double energy_lRatio = -999;
	int shower_axis_delta = -999;
	double vertex_radius = -999;
	int max_rms_shower_width = -999;
	double layer_min_energy = -999;
	int track_X_clusters = -999;
	int track_Y_clusters = -999;
	int track_missingHit_X = -999;
	int track_missingHit_Y = -999;
	int STK_BGO_delta_track = -999;
	int STK_BGO_delta_position = -999;
	int STK_PSD_delta_position = -999;
	double PSD_bar_min_energy_release = -999;
	double PSD_charge = -999;
	double PSD_charge_sum = -999;
	double STK_charge = -999;
};

struct active_cuts
{
	bool nBarLayer13 = false;
	bool maxRms = false;
	bool track_selection = false;
	bool psd_stk_match = false;
	bool psd_charge = false;
	bool stk_charge = false;
	unsigned int nActiveCuts = 0;
};

struct logger_cuts
{
	double sum_rms_min = -999;
	double sum_rms_max = -999;
};

struct logger_active_cuts
{
	bool sum_rms = false;
	bool trigger_only = false;
	bool bgo_only = false;
	bool all_cuts = false;
	unsigned int nActiveCuts = 0;
};

class config
{
public:
	config(
		const std::string working_dir,
		const bool mc);
	~config(){};
	std::vector<float> GetEnergyBinning();
	void PrintActiveFilters();
	const double GetMinEnergyRange();
	const double GetMaxEnergyRange();
	const double GetBGOLayerMinEnergy();
	const double GetPSDBarMinEnergy();
	const cuts_conf GetCutsConfigValues();
	const logger_cuts GetLoggerCutsConfigValues();
	const active_cuts GetActiveCuts();
	const logger_active_cuts GetLoggerActiveCuts();

private:
	std::string parse_config_file(
		std::string wd,
		std::string config_file);
	void get_config_info(std::string parsed_config);
	void get_logger_config_info(std::string parsed_config);
	cuts_conf cuts;
	active_cuts a_cuts;
	logger_cuts log_cuts;
	logger_active_cuts log_a_cuts;
	std::size_t n_bins;
	std::vector<float> energy_binning;
};

#endif