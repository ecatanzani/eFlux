#include "config.h"

config::config(
	const std::string working_dir,
	const bool mc)
{
	std::string config_file_name;
	mc ? config_file_name = "mc_config.conf" : config_file_name = "data_config.conf";
	get_config_info(parse_config_file(working_dir, config_file_name));
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
		// Load cuts variables
		if (!strcmp(tmp_str.c_str(), "min_event_energy"))
		{
			input_stream >> tmp_str;
			cuts.min_event_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "max_event_energy"))
		{
			input_stream >> tmp_str;
			cuts.max_event_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "energy_lRatio"))
		{
			input_stream >> tmp_str;
			cuts.energy_lRatio = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "shower_axis_delta"))
		{
			input_stream >> tmp_str;
			cuts.shower_axis_delta = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "max_rms_shower_width"))
		{
			input_stream >> tmp_str;
			cuts.max_rms_shower_width = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "layer_min_energy"))
		{
			input_stream >> tmp_str;
			cuts.layer_min_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "generation_vertex_radius"))
		{
			input_stream >> tmp_str;
			cuts.vertex_radius = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "track_X_clusters"))
		{
			input_stream >> tmp_str;
			cuts.track_X_clusters = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "track_Y_clusters"))
		{
			input_stream >> tmp_str;
			cuts.track_Y_clusters = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "track_missingHit_X"))
		{
			input_stream >> tmp_str;
			cuts.track_missingHit_X = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "track_missingHit_Y"))
		{
			input_stream >> tmp_str;
			cuts.track_missingHit_Y = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "STK_BGO_delta_track"))
		{
			input_stream >> tmp_str;
			cuts.STK_BGO_delta_track = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "STK_BGO_delta_position"))
		{
			input_stream >> tmp_str;
			cuts.STK_BGO_delta_position = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "STK_PSD_delta_position"))
		{
			input_stream >> tmp_str;
			cuts.STK_PSD_delta_position = stoi(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "PSD_bar_min_energy_release"))
		{
			input_stream >> tmp_str;
			cuts.PSD_bar_min_energy_release = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "PSD_charge"))
		{
			input_stream >> tmp_str;
			cuts.PSD_charge = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "PSD_charge_sum"))
		{
			input_stream >> tmp_str;
			cuts.PSD_charge_sum = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "PSD_charge_no_match"))
		{
			input_stream >> tmp_str;
			cuts.PSD_charge_no_match = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "STK_charge_upper"))
		{
			input_stream >> tmp_str;
			cuts.STK_charge_upper = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "STK_charge_medium"))
		{
			input_stream >> tmp_str;
			cuts.STK_charge_medium = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "STK_charge_lower"))
		{
			input_stream >> tmp_str;
			cuts.STK_charge_lower = stod(tmp_str, &sz);
		}


		// Load cuts
		if (!strcmp(tmp_str.c_str(), "nBarLayer13"))
		{
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
			{
				a_cuts.nBarLayer13 = true;
				++a_cuts.nActiveCuts;
			}
		}
		if (!strcmp(tmp_str.c_str(), "maxRms"))
		{
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
			{
				a_cuts.maxRms = true;
				++a_cuts.nActiveCuts;
			}
		}
		if (!strcmp(tmp_str.c_str(), "track_selection"))
		{
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
			{
				a_cuts.track_selection = true;
				++a_cuts.nActiveCuts;
			}
		}
		if (!strcmp(tmp_str.c_str(), "psd_stk_match"))
		{
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
			{
				a_cuts.psd_stk_match = true;
				++a_cuts.nActiveCuts;
			}
		}
		if (!strcmp(tmp_str.c_str(), "psd_charge"))
		{
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
			{
				a_cuts.psd_charge = true;
				++a_cuts.nActiveCuts;
			}
		}
		if (!strcmp(tmp_str.c_str(), "stk_charge"))
		{
			input_stream >> tmp_str;
			if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
			{
				a_cuts.stk_charge = true;
				++a_cuts.nActiveCuts;
			}
		}
	}
}

void config::PrintActiveFilters()
{
	std::string status_BarLayer13 = a_cuts.nBarLayer13 ? "ON" : "OFF";
	std::string status_maxRms = a_cuts.maxRms ? "ON" : "OFF";
	std::string status_track_selection = a_cuts.track_selection ? "ON" : "OFF";
	std::string status_psd_stk_match = a_cuts.psd_stk_match ? "ON" : "OFF";
	std::string status_psd_charge = a_cuts.psd_charge ? "ON" : "OFF";
	std::string status_stk_charge = a_cuts.stk_charge ? "ON" : "OFF";

	std::cout << "\n\n**** Filter Status ****\n";
	std::cout << "***********************\n\n";
	std::cout << "nBarLayer13 cut: " << status_BarLayer13 << std::endl;
	std::cout << "maxRms cut: " << status_maxRms << std::endl;
	std::cout << "track selection cut: " << status_track_selection << std::endl;
	std::cout << "PSD-STK match cut: " << status_psd_stk_match << std::endl;
	std::cout << "PSD charge cut: " << status_psd_charge << std::endl;
	std::cout << "STK charge cut: " << status_stk_charge << std::endl;
	std::cout << "\n***********************\n\n";
}

const double config::GetMinEnergyRange()
{
	return cuts.min_event_energy;
}

const double config::GetMaxEnergyRange()
{
	return cuts.max_event_energy;
}

const double config::GetBGOLayerMinEnergy()
{
	return cuts.layer_min_energy;
}

const double config::GetPSDBarMinEnergy()
{
	return cuts.PSD_bar_min_energy_release;
}

const cuts_conf config::GetCutsConfigValues()
{
	return cuts;
}

const active_cuts config::GetActiveCuts()
{
	return a_cuts;
}
