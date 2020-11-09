#ifndef HISTOS_H
#define HISTOS_H

#include <vector>
#include <memory>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

class histos
{
public:
	histos(std::vector<float> energy_bins);
	~histos(){};
	const int GetEnergyBin(const double energy);
	void FillTrigger(
		const double energy,
		const double energy_w);
	void FillBGO(
		const double raw_energy,
		const double corr_energy,
		const double simu_energy_w,
		const int reco_bidx,
		const std::vector<double> *fracLayer,
		const std::vector<double> *eLayer,
		const std::vector<double> *rmsLayer,
		const std::vector<double> *energy_1R_radius,
		const std::vector<double> *energy_2R_radius,
		const std::vector<double> *energy_3R_radius,
		const std::vector<double> *energy_5R_radius,
		const double sumRms,
		const double fracLast,
		const double fracLast_13,
		const int lastBGOLayer,
		const int nBGOentries,
		const double BGOrec_slopeX,
		const double BGOrec_slopeY,
		const double BGOrec_interceptX,
		const double BGOrec_interceptY,
		const double cosine_STK);
	void FillSumRmsCosine(
		const double energy,
		const double energy_w,
		const int reco_bidx,
		const double sum_rms,
		const double cosine);
	void FillSumRmsFLast(
		const double energy,
		const double energy_w,
		const int reco_bidx,
		const double sum_rms,
		const double last_energy_fraction,
		const double energy_fraction_13);
	void FillGeometric(
		const double energy,
		const double energy_w);
	void FillGeometricMaxLayer(
		const double energy,
		const double energy_w);
	void FillGeometricMaxBar(
		const double energy,
		const double energy_w);
	void FillGeometricBGOTrack(
		const double energy,
		const double energy_w);
	void FillGeometricBGOFiducial(
		const double energy,
		const double energy_w);
	void FillGeometricAll(
		const double energy,
		const double energy_w);
	void FillMaxLayer(
		const double energy,
		const double energy_w);
	void FillMaxBar(
		const double energy,
		const double energy_w);
	void FillBGOTrack(
		const double energy,
		const double energy_w);
	void FillBGOFiducial(
		const double energy,
		const double energy_w);
	void FillBGOFiducialBarL13(
		const double energy,
		const double energy_w);
	void FillBGOFiducialMaxRms(
		const double energy,
		const double energy_w);
	void FillBGOFiducialTrack(
		const double energy,
		const double energy_w);
	void FillBGOFiducialPsdStk(
		const double energy,
		const double energy_w);
	void FillBGOFiducialPsdCharge(
		const double energy,
		const double energy_w);
	void FillBGOFiducialStkCharge(
		const double energy,
		const double energy_w);
	void FillBGOFiducialAll(
		const double energy,
		const double energy_w);
	void FillNUD(
		const double energy_w,
		const std::vector<double> *nud_adc,
		const double nud_total_adc,
		const double nud_max_adc,
		const double nud_max_channel_id);
	void FillPsdCharge(
		const double energy_w,
		const double psd_charge_x,
		const double psd_charge_y,
		const double psd_charge,
		const bool selected=false);
	void FillStkCharge(
		const double energy_w,
		const double stk_charge_x,
		const double stk_charge_y,
		const double stk_charge,
		const bool selected=false);
	void FillAllCut(
		const double energy,
		const double energy_w);
	void FillClassifier(
		const double energy,
		const double energy_w,
		const double xtrl);
	void WriteCore(TFile* outfile);

protected:
	
	std::vector<float> logEBins;
	std::vector<float> cosine_bins;
	std::vector<float> sumRms_bins;
	std::vector<float> xtrl_bins;
	std::vector<float> flast_binning;
	bool simu = false;

	// Init
	void init_trigger_histos();
	void init_preselection_histos();
	void init_geometric_histos();
	void init_BGOfiducial_histos();
	void init_BGO_histos();
	void init_xtrl_histos();
	void init_ep_histos();
	void init_psd_charge_histos();
	void init_stk_charge_histos();
	void init_nud_histos();

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

	// Preselection - BGO fiducial cuts
	std::unique_ptr<TH1D> h_BGOfiducial_nBarLayer13_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_maxRms_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_track_selection_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_psd_stk_match_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_psd_charge_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_stk_charge_cut;
	std::unique_ptr<TH1D> h_BGOfiducial_all_cut;

	// BGO histos
	std::unique_ptr<TH1D> h_BGOrec_energy;
	std::unique_ptr<TH1D> h_BGOrec_corr_energy;
	std::unique_ptr<TH1D> h_BGOrec_layer_energy_diff;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_layer_max_energy_ratio;
	std::vector<std::vector<std::unique_ptr<TH1D>>> h_BGOrec_layer_energy_ratio;
	std::vector<std::vector<std::unique_ptr<TH1D>>> h_BGOrec_layer_rms;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_sumrms;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_sumrms_weighted;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_sumrms_cosine;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_fraclast;
	std::vector<std::unique_ptr<TH2D>> h_BGOrec_fraclast_cosine;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_frac13;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_last_layer;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_hits;
	std::vector<std::vector<std::unique_ptr<TH1D>>> h_BGOrec_energy_frac_1R;
	std::vector<std::vector<std::unique_ptr<TH1D>>> h_BGOrec_energy_frac_2R;
	std::vector<std::vector<std::unique_ptr<TH1D>>> h_BGOrec_energy_frac_3R;
	std::vector<std::vector<std::unique_ptr<TH1D>>> h_BGOrec_energy_frac_5R;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_slopeX;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_slopeY;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_interceptX;
	std::vector<std::unique_ptr<TH1D>> h_BGOrec_interceptY;
	std::vector<std::unique_ptr<TH2D>> h_BGOrec_topMap;
	std::vector<std::unique_ptr<TH2D>> h_BGOrec_bottomMap;
	std::vector<std::unique_ptr<TH2D>> h_BGOrec_shower_profile;
	std::vector<std::unique_ptr<TH2D>> h_BGOrec_shower_profile_cosine_upto_09;
	std::vector<std::unique_ptr<TH2D>> h_BGOrec_shower_profile_cone_from_09;

	std::vector<std::unique_ptr<TH2D>> sumRms_cosine;
	std::unique_ptr<TH2D> sumRms_cosine_20_100;
	std::unique_ptr<TH2D> sumRms_cosine_100_250;
	std::unique_ptr<TH2D> sumRms_cosine_250_500;
	std::unique_ptr<TH2D> sumRms_cosine_500_1000;
	std::unique_ptr<TH2D> sumRms_cosine_1000_3000;
	std::unique_ptr<TH2D> sumRms_cosine_3000_10000;
	std::unique_ptr<TH2D> sumRms_cosine_10000_20000;

	std::unique_ptr<TH1D> h_xtrl_energy_int;
	std::unique_ptr<TH2D> h_xtrl;
	std::vector<std::unique_ptr<TH1D>> h_xtrl_bin;

	std::vector<std::unique_ptr<TH2D>> e_discrimination_last;
	std::unique_ptr<TH2D> e_discrimination_last_20_100;
	std::unique_ptr<TH2D> e_discrimination_last_100_250;
	std::unique_ptr<TH2D> e_discrimination_last_250_500;
	std::unique_ptr<TH2D> e_discrimination_last_500_1000;
	std::unique_ptr<TH2D> e_discrimination_last_1000_3000;
	std::unique_ptr<TH2D> e_discrimination_last_3000_10000;
	std::unique_ptr<TH2D> e_discrimination_last_10000_20000;

	std::vector<std::unique_ptr<TH2D>> e_discrimination;
	std::unique_ptr<TH2D> e_discrimination_20_100;
	std::unique_ptr<TH2D> e_discrimination_100_250;
	std::unique_ptr<TH2D> e_discrimination_250_500;
	std::unique_ptr<TH2D> e_discrimination_500_1000;
	std::unique_ptr<TH2D> e_discrimination_1000_3000;
	std::unique_ptr<TH2D> e_discrimination_3000_10000;
	std::unique_ptr<TH2D> e_discrimination_10000_20000;

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

	// NUD
	std::vector<std::unique_ptr<TH1D>> h_NUD_adc;
	std::unique_ptr<TH1D> h_NUD_total_adc;
	std::unique_ptr<TH1D> h_NUD_max_adc;
	std::unique_ptr<TH1D> h_NUD_max_channel;
};

#endif