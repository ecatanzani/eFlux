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
		
		const std::string GetLEWeights();
		const std::string GetMEWeights();
		const std::string GetHEWeights();

	private:
		std::string parse_config_file(const char* config_file);
		void get_config_info(const std::string parsed_config);

		std::string le_weights;
		std::string me_weights;
		std::string he_weights;
};

#endif