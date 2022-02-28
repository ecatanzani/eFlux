#include "config.h"

config::config(const char* bdt_config_file) {
	get_config_info(parse_config_file(bdt_config_file));
}

std::string config::parse_config_file(const char* config_file) {
	std::ifstream input_file(config_file);
	if (!input_file.is_open())
	{
		std::cerr << "\nInput config file not found [" << config_file << "]\n\n";
		exit(100);
	}
	std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

void config::get_config_info(const std::string parsed_config) {
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);

	while (input_stream >> tmp_str)
	{
		// Load cuts variables
		if (!strcmp(tmp_str.c_str(), "low_energy_weights")) input_stream >> le_weights;
		if (!strcmp(tmp_str.c_str(), "medium_energy_weights")) input_stream >> me_weights;
		if (!strcmp(tmp_str.c_str(), "high_energy_weights")) input_stream >> he_weights;
	}
}

const std::string config::GetLEWeights() {
	return le_weights;
}

const std::string config::GetMEWeights() {
	return me_weights;
}

const std::string config::GetHEWeights() {
	return he_weights;
}

void config::PrintWeights() {

	std::cout << "\n\n**** Training Weights ****\n";
	std::cout << "***********************\n\n";
	std::cout << "Low Energy  (10 GeV - 100 GeV) : " << le_weights << std::endl;
	std::cout << "Mean energy (100 GeV - 1 TeV)  : " << me_weights << std::endl;
	std::cout << "High Energy (1 TeV - 10 TeV)   : " << he_weights << std::endl;
	std::cout << "\n***********************\n";
}