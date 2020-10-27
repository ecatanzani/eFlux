#ifndef HISTOS_H
#define HISTOS_H

#include <vector>
#include <memory>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

#include "DmpFilterContainer.h"

class histos
{
public:
	histos(std::vector<float> energy_bins);
	~histos(){};
	void Fill(
		const filter_output &output,
		const std::vector<double> &fracLayer,
		const std::vector<double> &bgoRec_slope,
		const std::vector<double> &bgoRec_intercept,
		const psd_charge &extracted_psd_charge,
		const stk_charge &extracted_stk_charge,
		const double sumRMS,
		const double costheta,
		const bgo_classifiers &classifier,
		const double lastFracLayer,
		const double frac_layer_13,
		const double energy,
		const double energy_w = 1,
		const double simu_energy = -999,
		const bool simu_evt = false);
	void Write(TFile &outfile);

protected:
	std::vector<float> logEBins;
	std::vector<float> cosine_bins;
	std::vector<float> sumRms_bins;
	std::vector<float> xtrl_bins;
	std::vector<float> flast_binning;
	bool simu = false;

	// Init
	void init_preselection_histos();
	void init_geometric_histos();
	void init_BGOfiducial_histos();
	void init_BGOlayer_histos();
	void init_sumRms_cosine_histos();
	void init_xtrl_histos();
	void init_ep_histos();
	void init_psd_charge_histos();
	void init_stk_charge_histos();
	void init_time_histos();
	void evaluate_energy_ratio(
		const std::vector<double> fracLayer,
		const double energy_w = 1);
	void evaluate_top_bottom_position(
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept,
		const double energy_w = 1);
	void fill_sumRms_cosine_histo(
		const double sumRMS,
		const double costheta,
		const double energy,
		const double energy_w = 1);
	void fill_XTRL_histo(
		const double xtrl,
		const double energy,
		const double energy_w = 1);
	void fill_ep_histo(
		const double sumRMS,
		const double lastFracLayer,
		const double frac_layer_13,
		const double energy,
		const double energy_w = 1);

	std::unique_ptr<TH1D> h_trigger;

	// Preselection cuts
	std::unique_ptr<TH1D> h_geometric_cut;
	std::unique_ptr<TH1D> h_maxElayer_cut;
	std::unique_ptr<TH1D> h_maxBarLayer_cut;
	std::unique_ptr<TH1D> h_BGOTrackContainment_cut;
	std::unique_ptr<TH1D> h_BGO_fiducial_cut;
	std::unique_ptr<TH1D> h_all_cut;

	// Preselection - geometric cuts
	std::unique_ptr<TH1D> h_geometric_maxElayer_cut;
	std::unique_ptr<TH1D> h_geometric_maxBarLayer_cut;
	std::unique_ptr<TH1D> h_geometric_BGOTrackContainment_cut;
	std::unique_ptr<TH1D> h_geometric_BGO_fiducial_cut;
	std::unique_ptr<TH1D> h_geometric_all_cut;
	std::unique_ptr<TH1D> h_geo_BGOrec_slopeX;
	std::unique_ptr<TH1D> h_geo_BGOrec_slopeY;
	std::unique_ptr<TH1D> h_geo_BGOrec_interceptX;
	std::unique_ptr<TH1D> h_geo_BGOrec_interceptY;
	std::unique_ptr<TH2D> h_geo_BGOreco_topMap;
	std::unique_ptr<TH2D> h_geo_BGOreco_bottomMap;

	// Preselection - BGO fiducial cuts
	std::unique_ptr<TH1D> h_BGOfiducial_nBarLayer13_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_maxRms_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_track_selection_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_psd_stk_match_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_psd_charge_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_stk_charge_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_all_cut;

	// Energy histos
	std::unique_ptr<TH1D> h_BGOrec_energy;
	std::unique_ptr<TH1D> h_layer_max_energy_ratio;
	std::vector<std::unique_ptr<TH1D>> h_layer_energy_ratio;

	// sumRms - cosine correlation
	std::vector<std::unique_ptr<TH2D>> sumRms_cosine;
	std::unique_ptr<TH2D> sumRms_cosine_20_100;
	std::unique_ptr<TH2D> sumRms_cosine_100_250;
	std::unique_ptr<TH2D> sumRms_cosine_250_500;
	std::unique_ptr<TH2D> sumRms_cosine_500_1000;
	std::unique_ptr<TH2D> sumRms_cosine_1000_3000;
	std::unique_ptr<TH2D> sumRms_cosine_3000_10000;

	// xtrl
	std::unique_ptr<TH1D> h_xtrl_energy_int;
	std::unique_ptr<TH2D> h_xtrl;
	std::vector<std::unique_ptr<TH1D>> h_xtrl_bin;

	// particle ID
	std::unique_ptr<TH2D> e_discrimination_last;
	std::unique_ptr<TH2D> e_discrimination_last_20_100;
	std::unique_ptr<TH2D> e_discrimination_last_100_250;
	std::unique_ptr<TH2D> e_discrimination_last_250_500;
	std::unique_ptr<TH2D> e_discrimination_last_500_1000;
	std::unique_ptr<TH2D> e_discrimination_last_1000_3000;
	std::unique_ptr<TH2D> e_discrimination_last_3000_10000;

	std::unique_ptr<TH2D> e_discrimination;
	std::unique_ptr<TH2D> e_discrimination_20_100;
	std::unique_ptr<TH2D> e_discrimination_100_250;
	std::unique_ptr<TH2D> e_discrimination_250_500;
	std::unique_ptr<TH2D> e_discrimination_500_1000;
	std::unique_ptr<TH2D> e_discrimination_1000_3000;
	std::unique_ptr<TH2D> e_discrimination_3000_10000;

	// PSD charge
	std::unique_ptr<TH1D> h_psd_chargeX;
	std::unique_ptr<TH1D> h_psd_chargeY;
	std::unique_ptr<TH2D> h_psd_charge2D;
	std::unique_ptr<TH1D> h_psd_charge;
	std::unique_ptr<TH1D> h_psd_selected_chargeX;
	std::unique_ptr<TH1D> h_psd_selected_chargeY;
	std::unique_ptr<TH2D> h_psd_selected_charge2D;
	std::unique_ptr<TH1D> h_psd_selected_charge;

	// STK charge
	std::unique_ptr<TH1D> h_stk_chargeX;
	std::unique_ptr<TH1D> h_stk_chargeY;
	std::unique_ptr<TH2D> h_stk_charge2D;
	std::unique_ptr<TH1D> h_stk_charge;
	std::unique_ptr<TH1D> h_stk_selected_chargeX;
	std::unique_ptr<TH1D> h_stk_selected_chargeY;
	std::unique_ptr<TH2D> h_stk_selected_charge2D;
	std::unique_ptr<TH1D> h_stk_selected_charge;

	// Event time
	//std::unique_ptr<TH1D> h_second;
};

#endif