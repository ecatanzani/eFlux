#include "mc_histos.h"
#include "binning.h"
#include "DAMPE_geo_structure.h"

#include "TDirectory.h"

void mc_histos::init_geo_factor()
{
	h_geo_factor = std::make_unique<TH1D>(
		"h_geo_factor",
		"Energy Distribution of the geometric factor; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
}

void mc_histos::init_simu_incoming()
{
	h_incoming = std::make_unique<TH1D>(
		"h_incoming",
		"Energy Distribution of the incoming particles; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
}

void mc_histos::init_simu_energy_histos()
{
	ratio_line_binning = createLinearBinning(0, 1, 100);
	h_simu_energy = std::make_unique<TH1D>(
		"h_simu_energy",
		"Simu Energy; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_simu_energy_w = std::make_unique<TH1D>(
		"h_simu_energy_w",
		"Simu Energy Weight",
		100, 0, 100);
	h_energy_diff = std::make_unique<TH1D>(
		"h_energy_diff",
		"Simu vs Raw Reco BGO energy: Simu Energy - Raw Energy (GeV); counts",
		1000, 0, 100);
	h_energy_diff_corr = std::make_unique<TH1D>(
		"h_energy_diff_corr",
		"Simu vs Corrected Reco BGO energy: Simu Energy - Corrected Energy (GeV); counts",
		1000, -100, 100);
	h_energy_diff2D = std::make_unique<TH2D>(
		"h_energy_diff2D",
		"Energy Ratio; Real Energy (GeV); (Raw - Simu)/Simu",
		logEBins.size() - 1,
		&(logEBins[0]),
		ratio_line_binning.size() - 1,
		&(ratio_line_binning[0]));
	h_energy_diff2D_corr = std::make_unique<TH2D>(
		"h_energy_diff2D_corr",
		"Energy Ratio; Real Energy (GeV); (Corr - Simu)/Simu",
		logEBins.size() - 1,
		&(logEBins[0]),
		ratio_line_binning.size() - 1,
		&(ratio_line_binning[0]));
	h_energy_unfold = std::make_unique<TH2D>(
		"h_energy_unfold",
		"Energy Unfolding Matrix; Real Energy (GeV); Raw Energy (GeV)",
		logEBins.size() - 1,
		&(logEBins[0]),
		logEBins.size() - 1,
		&(logEBins[0]));
	h_energy_unfold_corr = std::make_unique<TH2D>(
		"h_energy_unfold_corr",
		"Energy Unfolding Matrix; Real Energy (GeV); Corr Energy (GeV)",
		logEBins.size() - 1,
		&(logEBins[0]),
		logEBins.size() - 1,
		&(logEBins[0]));
}

void mc_histos::init_simu_histos()
{
	h_simu_BGOrec_topX_vs_realX = std::make_unique<TH1D>(
		"h_simu_BGOrec_topX_vs_realX",
		"Simu X - BGOrec TOP X",
		100, -100, 100);
	h_simu_BGOrec_topY_vs_realY = std::make_unique<TH1D>(
		"h_simu_BGOrec_topY_vs_realY",
		"Simu Y - BGOrec TOP Y",
		100, -100, 100);
	h_simu_slopeX = std::make_unique<TH1D>(
		"h_simu_slopeX",
		"Simu Slope X",
		1000, -90, 90);
	h_simu_slopeY = std::make_unique<TH1D>(
		"h_simu_slopeY",
		"Simu Slope Y",
		1000, -90, 90);
	h_simu_interceptX = std::make_unique<TH1D>(
		"h_simu_interceptX",
		"Simu Intercept X",
		500, -500, 500);
	h_simu_interceptY = std::make_unique<TH1D>(
		"h_simu_interceptY",
		"Simu Intercept Y",
		500, -500, 500);
	h_simu_topMap = std::make_unique<TH2D>(
		"h_simu_topMap",
		"Simu BGO TOP Map",
		500, -500, 500,
		500, -500, 500);
	h_simu_bottomMap = std::make_unique<TH2D>(
		"h_simu_bottomMap",
		"Simu BGO BOTTOM Map",
		500, -500, 500,
		500, -500, 500);
}

void mc_histos::FillGeoFactor(
	const double energy,
	const double energy_w)
{
	h_geo_factor->Fill(energy, energy_w);
}

void mc_histos::FillIncoming(
	const double energy,
	const double energy_w)
{
	h_incoming->Fill(energy, energy_w);
}

void mc_histos::FillSimu(
	const double simu_energy,
	const double simu_energy_w,
	const double raw_energy,
	const double corr_energy,
	const double simuSlopeX,
	const double simuSlopeY,
	const double simuInterceptX,
	const double simuInterceptY,
	const double BGOrec_slopeX,
	const double BGOrec_slopeY,
	const double BGOrec_interceptX,
	const double BGOrec_interceptY)
{
	const double _GeV = 0.001;
	h_simu_energy->Fill(simu_energy * _GeV, simu_energy_w);
	h_simu_energy_w->Fill(simu_energy_w, simu_energy_w);
	h_energy_diff->Fill((simu_energy - raw_energy)*_GeV, simu_energy_w);
	h_energy_diff_corr->Fill((simu_energy - corr_energy)*_GeV, simu_energy_w);
	h_energy_diff2D->Fill(simu_energy * _GeV, (raw_energy -simu_energy) / simu_energy, simu_energy_w);
	h_energy_diff2D_corr->Fill(simu_energy * _GeV, (corr_energy - simu_energy) / simu_energy, simu_energy_w);
	h_energy_unfold->Fill(simu_energy * _GeV, raw_energy * _GeV, simu_energy_w);
	h_energy_unfold_corr->Fill(simu_energy * _GeV, corr_energy * _GeV, simu_energy_w);
	h_simu_slopeX->Fill(simuSlopeX, simu_energy_w);
	h_simu_slopeY->Fill(simuSlopeY, simu_energy_w);
	h_simu_interceptX->Fill(simuInterceptX, simu_energy_w);
	h_simu_interceptY->Fill(simuInterceptY, simu_energy_w);

	double simu_topX = simuSlopeX * BGO_TopZ + simuInterceptX;
	double simu_topY = simuSlopeY * BGO_TopZ + simuInterceptY;
	double simu_bottomX = simuSlopeX * BGO_BottomZ + simuInterceptX;
	double simu_bottomY = simuSlopeY * BGO_BottomZ + simuInterceptY;
	double BGOrec_topX = BGOrec_slopeX * BGO_TopZ + BGOrec_interceptX;
	double BGOrec_topY = BGOrec_slopeY * BGO_TopZ + BGOrec_interceptY;
	//double BGOrec_bottomX = BGOrec_slopeX * BGO_BottomZ + BGOrec_interceptX;
	//double BGOrec_bottomY = BGOrec_slopeY * BGO_BottomZ + BGOrec_interceptY;
	double spread_topX = simu_topX - BGOrec_topX;
	double spread_topY = simu_topY - BGOrec_topY;

	h_simu_BGOrec_topX_vs_realX->Fill(spread_topX, simu_energy_w);
	h_simu_BGOrec_topY_vs_realY->Fill(spread_topY, simu_energy_w);;
	h_simu_topMap->Fill(simu_topX, simu_topY, simu_energy_w);
	h_simu_bottomMap->Fill(simu_bottomX, simu_bottomY, simu_energy_w);
}

void mc_histos::Write(TFile* outfile)
{
	outfile->cd();

	auto simu_dir = outfile->mkdir("Simu");
	simu_dir->cd();
	
	h_incoming->Write();
	h_simu_energy->Write();
	h_simu_energy_w->Write();
	h_energy_diff->Write();
	h_energy_diff_corr->Write();
	h_energy_diff2D->Write();
	h_energy_diff2D_corr->Write();
	h_energy_unfold->Write();
	h_energy_unfold_corr->Write();
	h_simu_BGOrec_topX_vs_realX->Write();
	h_simu_BGOrec_topY_vs_realY->Write();
	h_simu_slopeX->Write();
	h_simu_slopeY->Write();
	h_simu_interceptX->Write();
	h_simu_interceptY->Write();
	h_simu_topMap->Write();
	h_simu_bottomMap->Write();
}