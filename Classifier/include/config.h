#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

struct train_vars
{
	bool all_vars = false;
	bool no_nud = false;
	bool nud_only = false;
};

struct evt_stat
{
	unsigned int signal_train_events = 0;
	unsigned int background_train_events = 0;
	unsigned int signal_test_events = 0;
	unsigned int background_test_events = 0;
};

class config
{
public:
	config(const std::string working_dir)
	{
		get_config_info(parse_config_file(working_dir, config_file_name));
	}
	~config(){};
	
	const unsigned int GetSignalTrainEvents();
	const unsigned int GetSignalTestEvents();
	const unsigned int GetBackgroundTrainEvents();
	const unsigned int GetBackgroundTestEvents();
	const train_vars GetVariableOptions();

private:
	std::string parse_config_file(
		std::string wd,
		std::string config_file);
	void get_config_info(std::string parsed_config);
	
	std::string config_file_name = "classifier.conf";
	train_vars vars;
	evt_stat events;

};

#endif