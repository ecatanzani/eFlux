#ifndef MC_HISTOS_H
#define MC_HISTOS_H

#include "histos.h"

class mc_histos : public histos
{
public:
	mc_histos(std::vector<float> energy_bins) : histos(energy_bins)
	{
		init_geo_factor();
		init_incoming();
		init_energy_histos();
		init_geo_histos();
	};
	~mc_histos(){};
	void FillMC(
		const filter_output &output,
		const std::vector<double> &fracLayer,
		const std::vector<double> &bgoRec_slope,
		const std::vector<double> &bgoRec_intercept,
		const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
		const psd_charge &extracted_psd_charge,
		const stk_charge &extracted_stk_charge,
		const double sumRMS,
		const double costheta,
		const bgo_classifiers &classifier,
		const double lastFracLayer,
		const double frac_layer_13,
		const double simu_energy,
		const double raw_energy,
		const double corr_energy,
		const double energy_w = 1);
	void WriteMC(TFile &outfile);

private:
	void init_geo_factor();
	void init_incoming();
	void init_energy_histos();
	void init_pregeo_histos();
	void init_geo_histos();
	void evaluate_top_bottom_position_mc(
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept,
		const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
		const double energy_w = 1);

	std::vector<float> ratio_line_binning;

	// Geometric factor histo
	std::unique_ptr<TH1D> h_geo_factor;

	// Incoming histo
	std::unique_ptr<TH1D> h_incoming;

	// Energy histos
	std::unique_ptr<TH1D> h_simu_energy;
	std::unique_ptr<TH1D> h_energy_diff;
	std::unique_ptr<TH2D> h_energy_diff2D;
	std::unique_ptr<TH2D> h_energy_unfold;

	// geo histos
	std::unique_ptr<TH1D> h_geo_BGOrec_topX_vs_realX;
	std::unique_ptr<TH1D> h_geo_BGOrec_topY_vs_realY;
	std::unique_ptr<TH1D> h_geo_real_slopeX;
	std::unique_ptr<TH1D> h_geo_real_slopeY;
	std::unique_ptr<TH1D> h_geo_real_interceptX;
	std::unique_ptr<TH1D> h_geo_real_interceptY;
	std::unique_ptr<TH2D> h_geo_real_topMap;
	std::unique_ptr<TH2D> h_geo_real_bottomMap;
};

#endif