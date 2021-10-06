#ifndef CONFIG_H
#define CONFIG_H

#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TChain.h"
struct train_vars
{
	bool all_vars = false;
	bool no_nud = false;
	bool nud_only = false;
};

struct evt_stat
{
	bool auto_train_events = false;
	unsigned int signal_train_events = 0;
	unsigned int background_train_events = 0;
	unsigned int signal_test_events = 0;
	unsigned int background_test_events = 0;
};

struct cuts {
	bool xtrl;
	double signal_min_xtrl_value {0};
	double signal_max_xtrl_value {0};
	double background_min_xtrl_value {0};
	double background_max_xtrl_value {0};
	std::string signal_string, background_string;
};

class config
{
public:
	config(
		const std::string working_dir,
		std::shared_ptr<TChain> signal_train,
		std::shared_ptr<TChain> background_train);
	~config(){};
	
	const unsigned int GetSignalTrainEvents();
	const unsigned int GetSignalTestEvents();
	const unsigned int GetBackgroundTrainEvents();
	const unsigned int GetBackgroundTestEvents();
	const std::string GetSignalTCuts();
	const std::string GetBackgroundTCuts();
	const bool GetXTRLCutStatus();
	const train_vars GetVariableOptions();
	void PrintVariableOptions();

private:
	std::string parse_config_file(
		std::string wd,
		std::string config_file);
	void get_config_info(std::string parsed_config);
	
	std::string config_file_name = "classifier.conf";
	train_vars vars;
	evt_stat events;
	cuts applied_cuts;

};

#endif