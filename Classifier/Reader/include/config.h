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

class config
{
public:
	config(const std::string working_dir, const std::string bdt_config_file);
	~config(){};

	void PrintActiveFilters();
	void PrintWeights();

	const double GetMinEnergyRange();
	const double GetMaxEnergyRange();
	const double GetBGOLayerMinEnergy();
	const double GetPSDBarMinEnergy();
	const cuts_conf GetCutsConfigValues();
	const active_cuts GetActiveCuts();

	/*
	const std::string GetLEWeights();
	const std::string GetMEWeights();
	const std::string GetHEWeights();
	
	const double GetLEClassifierCut();
	const double GetMEClassifierCut();
	const double GetHEClassifierCut();
	*/

	const std::string GetBDTWeights_10_100();
    const std::string GetBDTWeights_100_250();
    const std::string GetBDTWeights_250_500();
    const std::string GetBDTWeights_500_1000();
    const std::string GetBDTWeights_1000_3000();
    const std::string GetBDTWeights_3000();

	const double GetBDTCut_10_100();
    const double GetBDTCut_100_250();
    const double GetBDTCut_250_500();
    const double GetBDTCut_500_1000();
    const double GetBDTCut_1000_3000();
    const double GetBDTCut_3000();

private:
	std::string parse_config_file(std::string config_file);
	void get_config_info(std::string parsed_config);
	void get_local_config_info(std::string parsed_config);
	cuts_conf cuts;
	active_cuts a_cuts;

	/*
	std::string le_weights;
	std::string me_weights;
	std::string he_weights;
	*/

	std::string weights_10_100;
	std::string weights_100_250;
	std::string weights_250_500;
    std::string weights_500_1000;
    std::string weights_1000_3000;
    std::string weights_3000;

	/*
	double le_c_cut {0};
	double me_c_cut {0};
	double he_c_cut {0};
	*/

	double cut_10_100 {0};
    double cut_100_250 {0};
    double cut_250_500 {0};
    double cut_500_1000 {0};
    double cut_1000_3000 {0};
    double cut_3000 {0};
};

#endif