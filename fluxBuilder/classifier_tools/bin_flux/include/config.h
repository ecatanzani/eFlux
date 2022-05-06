#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>


class config {
	public:
		config(const char* bdt_config_file);
		~config(){};

		void PrintWeights();
		
		/*
		const std::string GetLEWeights();
		const std::string GetMEWeights();
		const std::string GetHEWeights();
		*/

		const std::string GetBDTWeights_10_100();
		const std::string GetBDTWeights_100_250();
		const std::string GetBDTWeights_250_500();
		const std::string GetBDTWeights_500_1000();
		const std::string GetBDTWeights_1000_3000();
		const std::string GetBDTWeights_3000();

	private:
		std::string parse_config_file(const char* config_file);
		void get_config_info(const std::string parsed_config);

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
};

#endif