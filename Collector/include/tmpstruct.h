#ifndef TMPSTRUCT_H
#define TMPSTRUCT_H

#include <tuple>
#include <vector>
#include <memory>

#include "anyoption.h"
#include "Dmp/DmpNudContainer.h"
#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpFilterContainer.h"
#include "Efficiency/efficiency.h"
#include "Preselection/preselection.h"

#include "TVector3.h"

struct _tmp_filter
{
	filter_output output;
	presel_output preselection_output;
	eff_output efficiency_output;
	best_track evt_best_track;
	psd_charge evt_psd_charge;
	stk_charge evt_stk_charge;
	bgo_classifiers evt_bgo_classifier;
	trigger_info evt_trigger_info;
};

struct _tmp_psd
{
	double SPD_STK_match_X_distance; 					
	double SPD_STK_match_Y_distance; 					
	double SPD_STK_match_X_distance_fiducial_volume;	
	double SPD_STK_match_Y_distance_fiducial_volume;
};

struct _tmp_stk
{
	std::vector<int> clusters_on_plane;
	std::vector<double> stkEcore1Rm;
	std::vector<unsigned int> nStkClu1Rm;
};
struct _tmp_bgo
{
	std::vector<double> layer_energies;
	std::vector<std::vector<double>> layer_bar_energies;
	std::vector<double> slope;
	std::vector<double> intercept;
	TVector3 trajectory2D;
	double sumrms;
	std::vector<double> sumrms_layer;
	std::vector<double> energy_fraction_layer;
	double energy_fraction_last_layer;
	double energy_fraction_13th_layer;
	int last_energy_layer;
	int hits;
	double rvalue;										
	double lvalue;					
	double maximum_shower_position;				
	double maximum_shower_position_norm;
	std::vector<double> t_bgo;
	std::vector<double> t_bgo_norm;
};

struct _tmp_nud
{
	std::vector<int> adc;
	int total_adc;
	int max_adc;
	int max_channel_ID;
};

extern std::shared_ptr<_tmp_filter> fillFilterTmpStruct(DmpFilterContainer &filter, efficiency &eff_filter, preselection &presel_filter);
extern std::shared_ptr<_tmp_psd> fillPSDTmpStruct(DmpFilterContainer &filter);
extern std::shared_ptr<_tmp_stk> fillSTKTmpStruct(DmpStkContainer &stkVault);
extern std::shared_ptr<_tmp_bgo> fillBGOTmpStruct(DmpBgoContainer &bgoVault);
extern std::shared_ptr<_tmp_nud> fillNUDTmpStruct(DmpNudContainer &nudVault);

#endif