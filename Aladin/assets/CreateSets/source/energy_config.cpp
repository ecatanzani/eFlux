#include "energy_config.h"

energy_config::energy_config(const std::string working_dir) {
    get_config_info(parse_config_file(working_dir, config_file_name));
} 

std::string energy_config::parse_config_file(const std::string wd, const std::string config_file) {
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

void energy_config::get_config_info(const std::string parsed_config) {
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);
	std::string::size_type sz;

	while (input_stream >> tmp_str) {
		if (!strcmp(tmp_str.c_str(), "COMMON_ENERGY_RANGE")) {
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "ON") || !strcmp(tmp_str.c_str(), "on"))
				common_energy_range = true;
			else if (!strcmp(tmp_str.c_str(), "OFF") || !strcmp(tmp_str.c_str(), "off"))
				common_energy_range = false;
			else {
				std::cerr << "\n\nError parsing energy config file\n\n";
				exit(100);
			}
		}
		if (common_energy_range) {
			if (!strcmp(tmp_str.c_str(), "min_event_energy")) {
				input_stream >> tmp_str;
				min_event_energy = stod(tmp_str, &sz);
			}
			if (!strcmp(tmp_str.c_str(), "max_event_energy")) {
				input_stream >> tmp_str;
				max_event_energy = stod(tmp_str, &sz);
			}
		}
		else {
			if (!strcmp(tmp_str.c_str(), "train_min_event_energy")) {
				input_stream >> tmp_str;
				train_min_event_energy = stod(tmp_str, &sz);
			}
			if (!strcmp(tmp_str.c_str(), "train_max_event_energy")) {
				input_stream >> tmp_str;
				train_max_event_energy = stod(tmp_str, &sz);
			}
			if (!strcmp(tmp_str.c_str(), "test_min_event_energy")) {
				input_stream >> tmp_str;
				test_min_event_energy = stod(tmp_str, &sz);
			}
			if (!strcmp(tmp_str.c_str(), "test_max_event_energy")) {
				input_stream >> tmp_str;
				test_max_event_energy = stod(tmp_str, &sz);
			}
		}
	}
}

const double energy_config::GetMinEvtEnergy() {
	if (common_energy_range)
		return min_event_energy;
	else 
		return std::min(train_min_event_energy, test_min_event_energy);
}
    
const double energy_config::GetMaxEvtEnergy() {
	if (common_energy_range)
		return max_event_energy;
	else
		return std::max(train_max_event_energy, test_max_event_energy);
}

const double energy_config::GetMinTrainEvtEnergy() {
	return train_min_event_energy;
}

const double energy_config::GetMaxTrainEvtEnergy() {
	return train_max_event_energy;
}

const double energy_config::GetMinTestEvtEnergy() {
	return test_min_event_energy;
}

const double energy_config::GetMaxTestEvtEnergy() {
	return test_max_event_energy;
}

const bool energy_config::IsEnergyRangeCommon() {
	return common_energy_range;
}