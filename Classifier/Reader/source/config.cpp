#include "config.h"

config::config(const std::string working_dir, const std::string bdt_config_file) {
	get_config_info(parse_config_file(working_dir+std::string("/data_config.conf")));
	get_local_config_info(parse_config_file(bdt_config_file));
}

std::string config::parse_config_file(std::string config_file)
{
	std::ifstream input_file(config_file.c_str());
	if (!input_file.is_open())
	{
		std::cerr << "\nInput config file not found [" << config_file << "]\n\n";
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
		if (!strcmp(tmp_str.c_str(), "STK_charge"))
		{
			input_stream >> tmp_str;
			cuts.STK_charge = stod(tmp_str, &sz);
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

void config::get_local_config_info(std::string parsed_config) {
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);
	std::string::size_type sz;

	/*
	while (input_stream >> tmp_str)
	{
		// Load cuts variables
		if (!strcmp(tmp_str.c_str(), "low_energy_weights")) input_stream >> le_weights;
		if (!strcmp(tmp_str.c_str(), "medium_energy_weights")) input_stream >> me_weights;
		if (!strcmp(tmp_str.c_str(), "high_energy_weights")) input_stream >> he_weights;
		if (!strcmp(tmp_str.c_str(), "low_energy_classifier_cut")) { 
			input_stream >> tmp_str;
			le_c_cut = stod(tmp_str, &sz); 
		}
		if (!strcmp(tmp_str.c_str(), "medium_energy_classifier_cut")) { 
			input_stream >> tmp_str;
			me_c_cut = stod(tmp_str, &sz); 
		}
		if (!strcmp(tmp_str.c_str(), "high_energy_classifier_cut")) { 
			input_stream >> tmp_str;
			he_c_cut = stod(tmp_str, &sz); 
		}
	}
	*/

	while (input_stream >> tmp_str)
	{
		if (!strcmp(tmp_str.c_str(), "weights_10_100"))             input_stream >> weights_10_100;
		if (!strcmp(tmp_str.c_str(), "weights_100_250"))            input_stream >> weights_100_250;
		if (!strcmp(tmp_str.c_str(), "weights_250_500"))            input_stream >> weights_250_500;
        if (!strcmp(tmp_str.c_str(), "weights_500_1000"))           input_stream >> weights_500_1000;
        if (!strcmp(tmp_str.c_str(), "weights_1000_3000"))          input_stream >> weights_1000_3000;
        if (!strcmp(tmp_str.c_str(), "weights_3000"))               input_stream >> weights_3000;

        if (!strcmp(tmp_str.c_str(), "cut_10_100"))         {input_stream >> tmp_str; cut_10_100 = stod(tmp_str, &sz);}
        if (!strcmp(tmp_str.c_str(), "cut_100_250"))        {input_stream >> tmp_str; cut_100_250 = stod(tmp_str, &sz);}
        if (!strcmp(tmp_str.c_str(), "cut_250_500"))        {input_stream >> tmp_str; cut_250_500 = stod(tmp_str, &sz);}
        if (!strcmp(tmp_str.c_str(), "cut_500_1000"))       {input_stream >> tmp_str; cut_500_1000 = stod(tmp_str, &sz);}
        if (!strcmp(tmp_str.c_str(), "cut_1000_3000"))      {input_stream >> tmp_str; cut_1000_3000 = stod(tmp_str, &sz);}
        if (!strcmp(tmp_str.c_str(), "cut_3000"))           {input_stream >> tmp_str; cut_3000 = stod(tmp_str, &sz);}   
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
	std::cout << "\n***********************\n";
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

/*
const std::string config::GetLEWeights() {
	return le_weights;
}

const std::string config::GetMEWeights() {
	return me_weights;
}

const std::string config::GetHEWeights() {
	return he_weights;
}

const double config::GetLEClassifierCut() {
	return le_c_cut;
}

const double config::GetMEClassifierCut() {
	return me_c_cut;
}

const double config::GetHEClassifierCut() {
	return he_c_cut;
}

*/

const std::string config::GetBDTWeights_10_100()
{
	return weights_10_100;
}

const std::string config::GetBDTWeights_100_250()
{
	return weights_100_250;
}

const std::string config::GetBDTWeights_250_500()
{
	return weights_250_500;
}

const std::string config::GetBDTWeights_500_1000()
{
	return weights_500_1000;
}

const std::string config::GetBDTWeights_1000_3000()
{
	return weights_1000_3000;
}

const std::string config::GetBDTWeights_3000()
{
	return weights_3000;
}

const double config::GetBDTCut_10_100()
{
	return cut_10_100;
}

const double config::GetBDTCut_100_250()
{
	return cut_100_250;
}

const double config::GetBDTCut_250_500()
{
	return cut_250_500;
}

const double config::GetBDTCut_500_1000()
{
	return cut_500_1000;
}

const double config::GetBDTCut_1000_3000()
{
	return cut_1000_3000;
}

const double config::GetBDTCut_3000()
{
	return cut_3000;
}

/*
void config::PrintWeights() {

	std::cout << "\n\n**** Training Weights ****\n";
	std::cout << "***********************\n\n";
	std::cout << "Low Energy  (10 GeV - 100 GeV) : " << le_weights << std::endl;
	std::cout << "Mean energy (100 GeV - 1 TeV)  : " << me_weights << std::endl;
	std::cout << "High Energy (1 TeV - 10 TeV)   : " << he_weights << std::endl;
	std::cout << "Low Energy classifier cut  : " << le_c_cut << std::endl;
	std::cout << "Mean energy classifier cut : " << me_c_cut << std::endl;
	std::cout << "High Energy classifier cut : " << he_c_cut << std::endl;
	std::cout << "\n***********************\n";
}
*/

void config::PrintWeights() {

	std::cout << "\n\n**** Training Weights ****\n";
	std::cout << "***********************\n\n";
	std::cout << "Weights 10 GeV - 100 GeV: " << weights_10_100 << std::endl;
	std::cout << "Weights 100 GeV - 250 GeV: " << weights_100_250 << std::endl;
	std::cout << "Weights 250 GeV - 500 GeV: " << weights_250_500 << std::endl;
	std::cout << "Weights 500 GeV - 1 TeV: " << weights_500_1000 << std::endl;
	std::cout << "Weights 1 TeV - 3 TeV: " << weights_1000_3000 << std::endl;
	std::cout << "Weights 3 TeV - 10 TeV: " << weights_3000 << std::endl;


	std::cout << "Energy classifier cut  10 GeV - 100 GeV: " << cut_10_100 << std::endl;
	std::cout << "Energy classifier cut 100 GeV - 250 GeV: " << cut_100_250 << std::endl;
	std::cout << "Energy classifier cut 250 GeV - 500 GeV: " << cut_250_500 << std::endl;
	std::cout << "Energy classifier cut 500 GeV - 1 TeV: " << cut_500_1000 << std::endl;
	std::cout << "Energy classifier cut 1 TeV - 3 TeV: " << cut_1000_3000 << std::endl;
	std::cout << "Energy classifier cut 3 TeV - 10 TeV: " << cut_3000 << std::endl;
	std::cout << "\n***********************\n";
}