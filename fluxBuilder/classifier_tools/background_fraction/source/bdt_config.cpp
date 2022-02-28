#include "bdt_config.h"

bdt_config::bdt_config(const char* config_file) {
    get_config_info(parse_config_file(config_file));
}

std::string bdt_config::parse_config_file(const char* config_file_path) {
	std::ifstream input_file(config_file_path);
	if (!input_file.is_open()) {
		std::cerr << "\nInput config file not found [" << config_file_path << "]\n\n";
		exit(100);
	}
	std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

void bdt_config::get_config_info(const std::string parsed_config) {
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);
	std::string::size_type sz;

	while (input_stream >> tmp_str) {
		if (!strcmp(tmp_str.c_str(), "low_energy_classifier_cut")) { 
			input_stream >> tmp_str;
			le_c_cut = stod(tmp_str, &sz); 
		}
		if (!strcmp(tmp_str.c_str(), "medium_energy_classifier_cut")) { 
			input_stream >> tmp_str;
			me_c_cut = stod(tmp_str, &sz); 
		}
		if (!strcmp(tmp_str.c_str(), "high_energy_classifier_cut")) { 
			input_stream >> tmp_str;
			he_c_cut = stod(tmp_str, &sz); 
		}
	}
}

const double bdt_config::GetLowEnergyBDTCut() {
    return le_c_cut;
}

const double bdt_config::GetMidEnergyBDTCut() {
    return me_c_cut;
}

const double bdt_config::GetHighEnergyBDTCut() {
    return he_c_cut;
}