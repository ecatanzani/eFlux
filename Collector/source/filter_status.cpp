#include "acceptance.h"

void print_filter_status(data_active_cuts active_cuts)
{

	std::string status_BarLayer13 = active_cuts.nBarLayer13 ? "ON" : "OFF";
	std::string status_maxRms = active_cuts.maxRms ? "ON" : "OFF";
	std::string status_track_selection = active_cuts.track_selection ? "ON" : "OFF";
	std::string status_psd_stk_match = active_cuts.psd_stk_match ? "ON" : "OFF";
	std::string status_psd_charge = active_cuts.psd_charge ? "ON" : "OFF";
	std::string status_stk_charge = active_cuts.stk_charge ? "ON" : "OFF";
	std::string status_xtrl = active_cuts.xtrl ? "ON" : "OFF";

	std::cout << "\n\n**** Filter Status ****\n";
	std::cout << "***********************\n\n";
	std::cout << "nBarLayer13 cut: " << status_BarLayer13 << std::endl;
	std::cout << "maxRms cut: " << status_maxRms << std::endl;
	std::cout << "track selection cut: " << status_track_selection << std::endl;
	std::cout << "PSD-STK match cut: " << status_psd_stk_match << std::endl;
	std::cout << "PSD charge cut: " << status_psd_charge << std::endl;
	std::cout << "STK charge cut: " << status_stk_charge << std::endl;
	std::cout << "xtrl cut: " << status_xtrl;
	std::cout << "\n\n***********************\n\n";
}