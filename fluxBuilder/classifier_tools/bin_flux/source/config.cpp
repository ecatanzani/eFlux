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

	/*
	while (input_stream >> tmp_str)
	{
		// Load cuts variables
		if (!strcmp(tmp_str.c_str(), "low_energy_weights")) input_stream >> le_weights;
		if (!strcmp(tmp_str.c_str(), "medium_energy_weights")) input_stream >> me_weights;
		if (!strcmp(tmp_str.c_str(), "high_energy_weights")) input_stream >> he_weights;
	}
	*/

	while (input_stream >> tmp_str)
	{
		if (!strcmp(tmp_str.c_str(), "weights_10_100"))             input_stream >> weights_10_100;
		if (!strcmp(tmp_str.c_str(), "weights_100_250"))            input_stream >> weights_100_250;
		if (!strcmp(tmp_str.c_str(), "weights_250_500"))            input_stream >> weights_250_500;
        if (!strcmp(tmp_str.c_str(), "weights_500_1000"))           input_stream >> weights_500_1000;
        if (!strcmp(tmp_str.c_str(), "weights_1000_3000"))          input_stream >> weights_1000_3000;
        if (!strcmp(tmp_str.c_str(), "weights_3000"))               input_stream >> weights_3000;
	}
}

/*
const std::string config::GetLEWeights() {
	return le_weights;
}

const std::string config::GetMEWeights() {
	return me_weights;
}

const std::string config::GetHEWeights() {
	return he_weights;
}
*/

const std::string config::GetBDTWeights_10_100()
{
	return weights_10_100;
}

const std::string config::GetBDTWeights_100_250()
{
	return weights_100_250;
}

const std::string config::GetBDTWeights_250_500()
{
	return weights_250_500;
}

const std::string config::GetBDTWeights_500_1000()
{
	return weights_500_1000;
}

const std::string config::GetBDTWeights_1000_3000()
{
	return weights_1000_3000;
}

const std::string config::GetBDTWeights_3000()
{
	return weights_3000;
}

/*
void config::PrintWeights() {

	std::cout << "\n\n**** Training Weights ****\n";
	std::cout << "***********************\n\n";
	std::cout << "Low Energy  (10 GeV - 100 GeV) : " << le_weights << std::endl;
	std::cout << "Mean energy (100 GeV - 1 TeV)  : " << me_weights << std::endl;
	std::cout << "High Energy (1 TeV - 10 TeV)   : " << he_weights << std::endl;
	std::cout << "\n***********************\n";
}
*/

void config::PrintWeights() {

	std::cout << "\n\n**** Training Weights ****\n";
	std::cout << "***********************\n\n";
	std::cout << "Weights 10 GeV - 100 GeV: " << weights_10_100 << std::endl;
	std::cout << "Weights 100 GeV - 250 GeV: " << weights_100_250 << std::endl;
	std::cout << "Weights 250 GeV - 500 GeV: " << weights_250_500 << std::endl;
	std::cout << "Weights 500 GeV - 1 TeV: " << weights_500_1000 << std::endl;
	std::cout << "Weights 1 TeV - 3 TeV: " << weights_1000_3000 << std::endl;
	std::cout << "Weights 3 TeV - 10 TeV: " << weights_3000 << std::endl;
	std::cout << "\n***********************\n";
}