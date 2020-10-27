#include "mc_histos.h"
#include "binning.h"

void mc_histos::init_geo_factor()
{
	h_geo_factor = std::make_unique<TH1D>(
		"h_geo_factor",
		"Energy Distribution of the geometric factor; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_geo_factor->Sumw2();
}

void mc_histos::init_incoming()
{
	h_incoming = std::make_unique<TH1D>(
		"h_incoming",
		"Energy Distribution of the incoming particles; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_incoming->Sumw2();
}

void mc_histos::init_energy_histos()
{
	ratio_line_binning = createLinearBinning(0, 1, 100);
	h_simu_energy = std::make_unique<TH1D>(
		"h_simu_energy",
		"Simu Energy; Real Energy (GeV); counts",
		logEBins.size() - 1,
		&(logEBins[0]));
	h_simu_energy->Sumw2();
	h_energy_diff = std::make_unique<TH1D>(
		"h_energy_diff",
		"Simu vs Corrected Reco BGO energy: Real Energy - Corrected Energy (GeV); counts",
		100, 0, 1);
	h_energy_diff->Sumw2();
	h_energy_diff2D = std::make_unique<TH2D>(
		"h_energy_diff2D",
		"Energy Ratio; Real Energy (GeV); (Real - Raw)/Raw",
		logEBins.size() - 1,
		&(logEBins[0]),
		ratio_line_binning.size() - 1,
		&(ratio_line_binning[0]));
	h_energy_diff2D->Sumw2();
	h_energy_unfold = std::make_unique<TH2D>(
		"h_energy_unfold",
		"Energy Unfolding Matrix; Real Energy (GeV); Raw Energy (GeV)",
		logEBins.size() - 1,
		&(logEBins[0]),
		logEBins.size() - 1,
		&(logEBins[0]));
	h_energy_unfold->Sumw2();
}

void mc_histos::init_geo_histos()
{
	h_geo_BGOrec_topX_vs_realX = std::make_unique<TH1D>(
		"h_geo_BGOrec_topX_vs_realX",
		"Real X - BGOrec TOP X",
		100, -100, 100);
	h_geo_BGOrec_topX_vs_realX->Sumw2();
	h_geo_BGOrec_topY_vs_realY = std::make_unique<TH1D>(
		"h_geo_BGOrec_topY_vs_realY",
		"Real Y - BGOrec TOP Y",
		100, -100, 100);
	h_geo_BGOrec_topY_vs_realY->Sumw2();
	h_geo_real_slopeX = std::make_unique<TH1D>(
		"h_geo_real_slopeX",
		"Real Slope X",
		1000, -90, 90);
	h_geo_real_slopeX->Sumw2();
	h_geo_real_slopeY = std::make_unique<TH1D>(
		"h_geo_real_slopeY",
		"Real Slope Y",
		1000, -90, 90);
	h_geo_real_slopeY->Sumw2();
	h_geo_real_interceptX = std::make_unique<TH1D>(
		"h_geo_real_interceptX",
		"Real Intercept X",
		500, -500, 500);
	h_geo_real_interceptX->Sumw2();
	h_geo_real_interceptY = std::make_unique<TH1D>(
		"h_geo_real_interceptY",
		"Real Intercept Y",
		500, -500, 500);
	h_geo_real_interceptY->Sumw2();
	h_geo_real_topMap = std::make_unique<TH2D>(
		"h_geo_real_topMap",
		"Real BGO TOP Map",
		500, -500, 500,
		500, -500, 500);
	h_geo_real_topMap->Sumw2();
	h_geo_real_bottomMap = std::make_unique<TH2D>(
		"h_geo_real_bottomMap",
		"Real BGO BOTTOM Map",
		500, -500, 500,
		500, -500, 500);
	h_geo_real_bottomMap->Sumw2();
}

void mc_histos::evaluate_top_bottom_position_mc(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
	const double energy_w)
{
	// Get the real simu position
	TVector3 orgPosition;
	orgPosition.SetX(simu_primaries->pv_x);
	orgPosition.SetY(simu_primaries->pv_y);
	orgPosition.SetZ(simu_primaries->pv_z);

	TVector3 orgMomentum;
	orgMomentum.SetX(simu_primaries->pvpart_px);
	orgMomentum.SetY(simu_primaries->pvpart_py);
	orgMomentum.SetZ(simu_primaries->pvpart_pz);

	std::vector<double> slope(2, 0);
	std::vector<double> intercept(2, 0);

	slope[0] = orgMomentum.Z() ? orgMomentum.X() / orgMomentum.Z() : -999;
	slope[1] = orgMomentum.Z() ? orgMomentum.Y() / orgMomentum.Z() : -999;
	intercept[0] = orgPosition.X() - slope[0] * orgPosition.Z();
	intercept[1] = orgPosition.Y() - slope[1] * orgPosition.Z();

	double real_topX = slope[0] * BGO_TopZ + intercept[0];
	double real_topY = slope[1] * BGO_TopZ + intercept[1];

	double real_bottomX = slope[0] * BGO_BottomZ + intercept[0];
	double real_bottomY = slope[1] * BGO_BottomZ + intercept[1];

	double reco_topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
	double reco_topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

	double reco_bottomX = bgoRec_slope[0] * BGO_BottomZ + bgoRec_intercept[0];
	double reco_bottomY = bgoRec_slope[1] * BGO_BottomZ + bgoRec_intercept[1];

	h_geo_real_slopeX->Fill(slope[0], energy_w);
	h_geo_real_slopeY->Fill(slope[1], energy_w);
	h_geo_BGOrec_slopeX->Fill(bgoRec_slope[0], energy_w);
	h_geo_BGOrec_slopeY->Fill(bgoRec_slope[1], energy_w);
	h_geo_real_interceptX->Fill(intercept[0], energy_w);
	h_geo_real_interceptY->Fill(intercept[1], energy_w);
	h_geo_BGOrec_interceptX->Fill(bgoRec_intercept[0], energy_w);
	h_geo_BGOrec_interceptY->Fill(bgoRec_intercept[1], energy_w);

	auto spread_topX = real_topX - reco_topX;
	auto spread_topY = real_topY - reco_topY;

	h_geo_BGOrec_topX_vs_realX->Fill(spread_topX, energy_w);
	h_geo_BGOrec_topY_vs_realY->Fill(spread_topY, energy_w);

	h_geo_real_topMap->Fill(real_topX, real_topY, energy_w);
	h_geo_BGOreco_topMap->Fill(real_bottomX, real_bottomY, energy_w);
	h_geo_real_bottomMap->Fill(reco_topX, reco_topY, energy_w);
	h_geo_BGOreco_bottomMap->Fill(reco_bottomX, reco_bottomY, energy_w);
}

void mc_histos::FillMC(
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
	const double energy_w)
{
	auto _GeV = 0.001;
	auto simu_energy_gev = simu_energy * _GeV;
	auto raw_energy_gev = raw_energy * _GeV;
	auto corr_energy_gev = corr_energy * _GeV;
	
	Fill(
		output,
		fracLayer,
		bgoRec_slope,
		bgoRec_intercept,
		extracted_psd_charge,
		extracted_stk_charge,
		sumRMS,
		costheta,
		classifier,
		lastFracLayer,
		frac_layer_13,
		corr_energy,
		energy_w,
		simu_energy,
		true);
	
	if (!output.out_energy_range)
	{
		if (output.geometric_before_trigger)
		{
			h_geo_factor->Fill(simu_energy_gev, energy_w);
			if (output.evt_triggered)
			{
				h_trigger->Fill(simu_energy_gev, energy_w);
				if (output.correct_bgo_reco)
				{
					h_incoming->Fill(simu_energy_gev, energy_w);
					evaluate_top_bottom_position_mc(
						bgoRec_slope,
						bgoRec_intercept,
						simu_primaries,
						energy_w);
					h_simu_energy->Fill(
						simu_energy_gev,
						energy_w);
					h_energy_diff->Fill(
						(simu_energy - corr_energy) / simu_energy,
						energy_w);
					h_energy_diff2D->Fill(
						simu_energy_gev,
						(simu_energy - corr_energy) / simu_energy,
						energy_w);
					h_energy_unfold->Fill(
						simu_energy_gev,
						corr_energy_gev,
						energy_w);
				}
			}
			else
				h_incoming->Fill(simu_energy_gev, energy_w);
		}
	}
}

void mc_histos::WriteMC(TFile &outfile)
{
	Write(outfile);
	h_geo_factor->Write();
	h_incoming->Write();

	outfile.cd("BGO_Rec");

	h_geo_BGOrec_topX_vs_realX->Write();
	h_geo_BGOrec_topY_vs_realY->Write();
	h_geo_real_slopeX->Write();
	h_geo_real_slopeY->Write();
	h_geo_real_interceptX->Write();
	h_geo_real_interceptY->Write();
	h_geo_real_topMap->Write();
	h_geo_real_bottomMap->Write();

	outfile.cd("BGO_Energy");

	h_simu_energy->Write();
	h_energy_diff->Write();
	h_energy_diff2D->Write();
	h_energy_unfold->Write();
}