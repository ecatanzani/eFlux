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
	config(const std::string working_dir);
	~config(){};

	const double GetMinEnergyRange();
	const double GetMaxEnergyRange();

private:
	std::string parse_config_file(std::string wd, std::string config_file);
	void get_config_info(std::string parsed_config);
	
	double min_evt_energy {0};
	double max_evt_energy {0};

};

#endif