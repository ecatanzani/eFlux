#include "config.h"

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

const train_vars config::GetVariableOptions()
{
	return vars;
}