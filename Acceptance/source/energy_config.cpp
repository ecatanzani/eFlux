#include "energy_config.h"

energy_config::energy_config(const std::string working_dir)
{
    const std::string local_path = "/Acceptance/config";
    const std::string config_file_name = "energy_config.conf";
    const auto tmp_config = working_dir.substr(0, working_dir.find("/Collector/config")) + local_path;
    get_config_info(parse_config_file(tmp_config, config_file_name));
	n_bins = energy_binning.size()-1;
	min_event_energy = energy_binning.front();
	max_event_energy = energy_binning.back();
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
	
	int introduction {4};
	int elm_in_row {5};
    std::vector<std::string> row;

	int column_counter {0};
    int line_counter {0};

	while (input_stream >> tmp_str) {

        if (!line_counter) {
			// This is the first line... we are not interested in it
			++column_counter;
			if (column_counter==introduction) {
				++line_counter;
				column_counter = 0;
			}
		}
		else {
			// This is a general line...
			row.push_back(tmp_str);
			++column_counter;
			
			if (column_counter == elm_in_row) {

				// The row of the binning has been completed... let's extract the info
				if (line_counter==1) 
					energy_binning.push_back(stod(row[2], &sz));
				energy_binning.push_back(stod(row.back(), &sz));

				// Reset
				column_counter = 0;
				++line_counter;
				row.clear();
			}
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
	std::cout << "\n****************************\n\n";
}

std::vector<double> energy_config::GetEnergyBinning()
{
	return energy_binning;
}

const double energy_config::GetMinEvtEnergy()
{
	return min_event_energy;
}

const double energy_config::GetMaxEvtEnergy()
{
	return max_event_energy;
}