#ifndef TMPSTRUCT_H
#define TMPSTRUCT_H

#include <vector>
#include <memory>

#include "anyoption.h"
#include "Dmp/DmpNudContainer.h"
#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpFilterContainer.h"

#include "TVector3.h"

struct _tmp_filter
{
	filter_output output;
	p_cuts preselection_output;
	best_track evt_best_track;
	psd_charge evt_psd_charge;
	stk_charge evt_stk_charge;
	bgo_classifiers evt_bgo_classifier;
	trigger_info evt_trigger_info;
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
};

struct _tmp_nud
{
	std::vector<int> adc;
	int total_adc;
	int max_adc;
	int max_channel_ID;
};

extern std::shared_ptr<_tmp_filter> fillFilterTmpStruct(DmpFilterContainer &filter);
extern std::shared_ptr<_tmp_bgo> fillBGOTmpStruct(DmpBgoContainer &bgoVault);
extern std::shared_ptr<_tmp_nud> fillNUDTmpStruct(DmpNudContainer &nudVault);

#endif