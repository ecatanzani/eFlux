#include "config.h"

#include <algorithm>

config::config(
	const std::string working_dir,
	std::shared_ptr<TChain> signal_train,
	std::shared_ptr<TChain> background_train)
{
	// Parse config file
	get_config_info(parse_config_file(working_dir, config_file_name));
	// Auto-set training events
	if (events.auto_train_events)
		events.signal_train_events = events.background_train_events = std::min(signal_train->GetEntries(applied_cuts.signal_string.c_str()), background_train->GetEntries(applied_cuts.background_string.c_str()));
}

std::string config::parse_config_file(
	std::string wd,
	std::string config_file)
{
	std::string configPath = wd + "/" + config_file;
	std::ifstream input_file(configPath.c_str());
	if (!input_file.is_open())
	{
		std::cerr << "\nInput config file not found [" << configPath << "]\n\n";
		exit(100);
	}
	std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

void config::get_config_info(std::string parsed_config)
{
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);
	std::string::size_type sz;

	while (input_stream >> tmp_str)
	{
		// Load vars config
		if (!strcmp(tmp_str.c_str(), "all_vars"))
		{
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
				vars.all_vars = true;
		}
		if (!strcmp(tmp_str.c_str(), "NUD_only"))
		{
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
				vars.nud_only = true;
		}
		if (!strcmp(tmp_str.c_str(), "no_NUD"))
		{
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
				vars.no_nud = true;
		}

		// Load train events
		if (!strcmp(tmp_str.c_str(), "auto_train_events"))
		{
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
				events.auto_train_events = true;
		}
		if (!events.auto_train_events)
		{
			if (!strcmp(tmp_str.c_str(), "signal_train_events"))
			{
				input_stream >> tmp_str;
				events.signal_train_events = stoul(tmp_str, &sz);
			}
			if (!strcmp(tmp_str.c_str(), "background_train_events"))
			{
				input_stream >> tmp_str;
				events.background_train_events = stoul(tmp_str, &sz);
			}
		}
		
		// Load test events
		if (!strcmp(tmp_str.c_str(), "signal_test_events"))
		{
			input_stream >> tmp_str;
			events.signal_test_events = stoul(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "background_test_events"))
		{
			input_stream >> tmp_str;
			events.background_test_events = stoul(tmp_str, &sz);
		}

		// Load cuts
		if (!strcmp(tmp_str.c_str(), "XTRL")) {
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "ON") || !strcmp(tmp_str.c_str(), "on"))
				applied_cuts.xtrl = true;
			else if (!strcmp(tmp_str.c_str(), "OFF") || !strcmp(tmp_str.c_str(), "off"))
				applied_cuts.xtrl = false;
			else {
				std::cerr << "\n\nError parsing energy config file\n\n";
				exit(100);
			}
		}
		if (applied_cuts.xtrl) {
			if (!strcmp(tmp_str.c_str(), "signal_min_value")) {
				input_stream >> tmp_str;
				applied_cuts.signal_min_xtrl_value = stod(tmp_str, &sz);
			}
			if (!strcmp(tmp_str.c_str(), "signal_max_value")) {
				input_stream >> tmp_str;
				applied_cuts.signal_max_xtrl_value = stod(tmp_str, &sz);
			}
			if (!strcmp(tmp_str.c_str(), "background_min_value")) {
				input_stream >> tmp_str;
				applied_cuts.background_min_xtrl_value = stod(tmp_str, &sz);
			}
			if (!strcmp(tmp_str.c_str(), "background_max_value")) {
				input_stream >> tmp_str;
				applied_cuts.background_max_xtrl_value = stod(tmp_str, &sz);
			}
			applied_cuts.signal_string = std::string("xtrl>") + std::to_string(applied_cuts.signal_min_xtrl_value) + std::string("&& xtrl<") + std::to_string(applied_cuts.signal_max_xtrl_value);
			applied_cuts.background_string = std::string("xtrl>") + std::to_string(applied_cuts.background_min_xtrl_value) + std::string("&& xtrl<") + std::to_string(applied_cuts.background_max_xtrl_value);
		}

		// fraclastlayer_norm removal
		if (applied_cuts.signal_string.empty()) {
			applied_cuts.signal_string = std::string("fraclastlayer_norm>=-100 && fraclastlayer_norm<=100");
			applied_cuts.background_string = std::string("fraclastlayer_norm>=-100 && fraclastlayer_norm<=100");
		}
		else {
			applied_cuts.signal_string += std::string(" && fraclastlayer_norm>=-100 && fraclastlayer_norm<=100");
			applied_cuts.background_string += std::string("&& fraclastlayer_norm>=-100 && fraclastlayer_norm<=100");
		}
	}
}

const unsigned int config::GetSignalTrainEvents()
{
	return events.signal_train_events;
}

const unsigned int config::GetSignalTestEvents()
{
	return events.signal_test_events;
}

const unsigned int config::GetBackgroundTrainEvents()
{
	return events.background_train_events;
}

const unsigned int config::GetBackgroundTestEvents()
{
	return events.background_test_events;
}

const bool config::GetXTRLCutStatus() {
	return applied_cuts.xtrl;
}

const std::string config::GetSignalTCuts() {
	return applied_cuts.signal_string;
}

const std::string config::GetBackgroundTCuts() {
	return applied_cuts.background_string;
}


const train_vars config::GetVariableOptions()
{
	return vars;
}

void config::PrintVariableOptions()
{
	std::string all_vars = vars.all_vars ? "ON" : "OFF";
	std::string no_nud = vars.no_nud ? "ON" : "OFF";
	std::string nud_only = vars.nud_only ? "ON" : "OFF";

	std::cout << "\n\n**** Variable Status ****\n";
	std::cout << "***********************\n\n";
	std::cout << "- all vars: " << all_vars << std::endl;
	std::cout << "- no_nud: " << no_nud << std::endl;
	std::cout << "- nud_only: " << nud_only << std::endl;
	std::cout << "\n***********************\n";

	if (events.auto_train_events)
	{
		std::cout << "\n\n Auto train events ACTIVATED - signal training events: " << events.signal_train_events;
		std::cout << "\n\n Auto train events ACTIVATED - background training events: " << events.background_train_events << std::endl << std::endl;
	}
}