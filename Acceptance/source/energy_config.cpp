#include "binning.h"
#include "energy_config.h"

energy_config::energy_config(const std::string working_dir)
{
    const std::string local_path = "/Acceptance/config";
    const std::string config_file_name = "energy_config.conf";
    const auto tmp_config = working_dir.substr(0, working_dir.find("/Collector/config")) + local_path;
    get_config_info(parse_config_file(tmp_config, config_file_name));
    energy_binning = createLogBinning(
		min_event_energy,
		max_event_energy,
		n_bins);
} 

std::string energy_config::parse_config_file(
	const std::string wd,
	const std::string config_file)
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

void energy_config::get_config_info(const std::string parsed_config)
{
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);
	std::string::size_type sz;

	while (input_stream >> tmp_str)
	{
		if (!strcmp(tmp_str.c_str(), "n_energy_bins"))
			input_stream >> n_bins;
		if (!strcmp(tmp_str.c_str(), "min_event_energy"))
		{
			input_stream >> tmp_str;
			min_event_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "max_event_energy"))
		{
			input_stream >> tmp_str;
			max_event_energy = stod(tmp_str, &sz);
		}
		// Energy filter
		if (!strcmp(tmp_str.c_str(), "min_event_energy_set"))
		{
			input_stream >> tmp_str;
			set_filter.min_event_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "max_event_energy_set"))
		{
			input_stream >> tmp_str;
			set_filter.max_event_energy = stod(tmp_str, &sz);
		}
	}
}

void energy_config::PrintActiveFilters()
{
	std::cout << "\n**** Energy Config File ****\n";
	std::cout << "****************************\n\n";
	std::cout << "Number of energy bins: " << n_bins << std::endl;
	std::cout << "Min event energy: " << min_event_energy << std::endl;
	std::cout << "Max event energy: " << max_event_energy << std::endl;
	std::cout << std::endl << "Set min event energy: " << set_filter.min_event_energy << std::endl;
	std::cout << "Set max event energy: " << set_filter.max_event_energy << std::endl;
	std::cout << "\n****************************\n\n";
}

std::vector<float> energy_config::GetEnergyBinning()
{
	return energy_binning;
}

const double energy_config::GetSetMinEvtEnergy()
{
	return set_filter.min_event_energy;
}

const double energy_config::GetSetMaxEvtEnergy()
{
	return set_filter.max_event_energy;
}