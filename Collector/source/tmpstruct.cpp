#include "tmpstruct.h"

std::shared_ptr<_tmp_filter> fillFilterTmpStruct(DmpFilterContainer &filter, efficiency &eff_filter)
{
	std::shared_ptr<_tmp_filter> _filter_res = std::make_shared<_tmp_filter>();	
	_filter_res->output = filter.GetFilterOutput();
	_filter_res->efficiency_preselection_output = eff_filter.GetEfficiencyOutput();
	_filter_res->preselection_output = filter.GetPreselectionOutput();
	_filter_res->evt_best_track = filter.GetBestTrack();
	_filter_res->evt_psd_charge = filter.GetPSDCharge();
	_filter_res->evt_stk_charge = filter.GetSTKCharge();
	_filter_res->evt_bgo_classifier = filter.GetClassifiers();
	_filter_res->evt_trigger_info = filter.GetTrigger();

	return _filter_res;
}

std::shared_ptr<_tmp_bgo> fillBGOTmpStruct(DmpBgoContainer &bgoVault)
{
	std::shared_ptr<_tmp_bgo> _bgo_res = std::make_shared<_tmp_bgo>();

	_bgo_res->layer_energies = bgoVault.GetLayerEnergies();
	_bgo_res->layer_bar_energies = bgoVault.GetLayerBarEnergies();
	_bgo_res->slope = bgoVault.GetBGOslope();
	_bgo_res->intercept = bgoVault.GetBGOintercept();
	_bgo_res->trajectory2D = bgoVault.GetBGOTrajectory2D();
	_bgo_res->sumrms = bgoVault.GetSumRMS();
	_bgo_res->sumrms_layer = bgoVault.GetRmsLayer();
	_bgo_res->energy_fraction_layer = bgoVault.GetFracLayer();
	_bgo_res->energy_fraction_last_layer = bgoVault.GetSingleFracLayer(bgoVault.GetLastEnergyLayer());
	_bgo_res->energy_fraction_13th_layer = bgoVault.GetSingleFracLayer(13);
	_bgo_res->last_energy_layer = bgoVault.GetLastEnergyLayer();
	_bgo_res->hits = bgoVault.GetNhits();
	return _bgo_res;
}

std::shared_ptr<_tmp_nud> fillNUDTmpStruct(DmpNudContainer &nudVault)
{
	std::shared_ptr<_tmp_nud> _nud_res = std::make_shared<_tmp_nud>();
	
	_nud_res->adc = nudVault.GetADC();
	_nud_res->total_adc = nudVault.GetTotalADC();
	_nud_res->max_adc = nudVault.GetMaxADC();
	_nud_res->max_channel_ID = nudVault.GetMaxChannelID();

	return _nud_res;
}