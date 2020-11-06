#ifndef MC_HISTOS_H
#define MC_HISTOS_H

#include "histos.h"

#include "TVector3.h"

class mc_histos : public histos
{
public:
	mc_histos(std::vector<float> energy_bins) : histos(energy_bins)
	{
		init_geo_factor();
		init_simu_incoming();
		init_simu_energy_histos();
		init_simu_histos();
	};
	~mc_histos(){};
	void FillGeoFactor(
		const double energy,
		const double energy_w);
	void FillIncoming(
		const double energy,
		const double energy_w);
	void FillSimu(
		const double simu_energy,
		const double simu_energy_w,
		const double corr_energy,
		const double simuSlopeX,
		const double simuSlopeY,
		const double simuInterceptX,
		const double simuInterceptY,
		const double BGOrec_slopeX,
		const double BGOrec_slopeY,
		const double BGOrec_interceptX,
		const double BGOrec_interceptY);
	void Write(TFile* outfile);
	
private:
	void init_geo_factor();
	void init_simu_incoming();
	void init_simu_energy_histos();
	void init_pregeo_histos();
	void init_simu_histos();

	std::vector<float> ratio_line_binning;

	// Geometric factor histo
	std::unique_ptr<TH1D> h_geo_factor;

	// Incoming histo
	std::unique_ptr<TH1D> h_incoming;

	// Energy histos
	std::unique_ptr<TH1D> h_simu_energy;
	std::unique_ptr<TH1D> h_simu_energy_w;
	std::unique_ptr<TH1D> h_energy_diff;
	std::unique_ptr<TH2D> h_energy_diff2D;
	std::unique_ptr<TH2D> h_energy_unfold;

	// Simu histos
	std::unique_ptr<TH1D> h_simu_BGOrec_topX_vs_realX;
	std::unique_ptr<TH1D> h_simu_BGOrec_topY_vs_realY;
	std::unique_ptr<TH1D> h_simu_slopeX;
	std::unique_ptr<TH1D> h_simu_slopeY;
	std::unique_ptr<TH1D> h_simu_interceptX;
	std::unique_ptr<TH1D> h_simu_interceptY;
	std::unique_ptr<TH2D> h_simu_topMap;
	std::unique_ptr<TH2D> h_simu_bottomMap;
};

#endif