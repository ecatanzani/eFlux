#include "config.h"

#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

config::config(const std::string energy_config_path) {
    get_config_info(parse_config_file(energy_config_path.substr(0, energy_config_path.find("/config")) + std::string("/Collector/assets/PreSelection/config"), std::string("cuts.conf")));
}

std::string config::parse_config_file(std::string wd,std::string config_file) {
	std::string configPath = wd + "/" + config_file;
	std::ifstream input_file(configPath.c_str());
	if (!input_file.is_open()){
		std::cerr << "\nInput config file not found [" << configPath << "]\n\n";
		exit(100);
	}
	std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

void config::get_config_info(std::string parsed_config)
{
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);
	std::string::size_type sz;

	while (input_stream >> tmp_str) {
		
        if (!strcmp(tmp_str.c_str(), "energy_lRatio")) {
			input_stream >> tmp_str;
			cuts.bgo_max_energy_ratio = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "shower_axis_delta")) {
			input_stream >> tmp_str;
			cuts.bgo_shower_axis_delta = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "max_rms_shower_width")) {
			input_stream >> tmp_str;
			cuts.bgo_shower_width = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "layer_min_energy")) {
			input_stream >> tmp_str;
			cuts.bgo_layer_min_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "track_X_clusters")) {
			input_stream >> tmp_str;
			cuts.track_X_clusters = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "track_Y_clusters")) {
			input_stream >> tmp_str;
			cuts.track_Y_clusters = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "track_missingHit_X")) {
			input_stream >> tmp_str;
			cuts.track_X_holes = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "track_missingHit_Y")) {
			input_stream >> tmp_str;
			cuts.track_Y_holes = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "STK_BGO_delta_track")) {
			input_stream >> tmp_str;
			cuts.STK_BGO_delta_track = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "STK_BGO_delta_position")) {
			input_stream >> tmp_str;
			cuts.STK_BGO_delta_position = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "STK_PSD_delta_position")) {
			input_stream >> tmp_str;
			cuts.STK_PSD_delta_position = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "PSD_bar_min_energy_release")) {
			input_stream >> tmp_str;
			cuts.psd_min_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "PSD_charge")) {
			input_stream >> tmp_str;
			cuts.PSD_single_charge = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "PSD_charge_sum")) {
			input_stream >> tmp_str;
			cuts.PSD_sharge_sum = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "STK_charge")) {
			input_stream >> tmp_str;
			cuts.STK_single_charge = stod(tmp_str, &sz);
        }
	}
}

const cuts_config config::GetCutsConfig() {
    return cuts;
}