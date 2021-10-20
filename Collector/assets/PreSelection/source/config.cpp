#include "config.h"

config::config(const std::string working_dir) {
	get_config_info(parse_config_file(working_dir, "preselection.conf"));
}

std::string config::parse_config_file(std::string wd, std::string config_file) {
	std::string configPath = wd + "/" + config_file;
	std::ifstream input_file(configPath.c_str());
	if (!input_file.is_open()) {
		std::cerr << "\nInput config file not found [" << configPath << "]\n\n";
		exit(100);
	}
	std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

void config::get_config_info(std::string parsed_config) {
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);
	std::string::size_type sz;

	while (input_stream >> tmp_str) {
		// Load cuts variables
		if (!strcmp(tmp_str.c_str(), "min_event_energy"))
		{
			input_stream >> tmp_str;
			min_evt_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "max_event_energy"))
		{
			input_stream >> tmp_str;
			max_evt_energy = stod(tmp_str, &sz);
		}
	}
}

const double config::GetMinEnergyRange() {
	return min_evt_energy;
}

const double config::GetMaxEnergyRange() {
	return max_evt_energy;
}