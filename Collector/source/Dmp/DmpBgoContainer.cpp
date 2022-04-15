#include "Dmp/DmpBgoContainer.h"

#include <cmath>
#include <iostream>

#include "TF1.h"
#include "TGraph.h"

void DmpBgoContainer::scanBGOHits(
	const std::shared_ptr<DmpEvtBgoHits> bgohits,
	const std::shared_ptr<DmpEvtBgoRec> bgorec,
	const double bgoTotalE,
	const double layerMinEnergy,
	const int nLayers)
{
	// Compute BGO shower axis
	if (slope[0] == -999 && slope[1] == -999)
	{
		slope[0] = bgorec->GetSlopeXZ();
		slope[1] = bgorec->GetSlopeYZ();
	}
	if (intercept[0] == -999 && intercept[1] == -999)
	{
		intercept[0] = bgorec->GetInterceptXZ();
		intercept[1] = bgorec->GetInterceptYZ();
	}
	// Get BGO reco 2D trajectory
	trajectoryDirection2D = bgorec->GetTrajectoryDirection2D();
	
	// Get the number of BGO hits
	nBgoHits = bgohits->GetHittedBarNumber();
	
	// Scan BGO hits
	for (int ihit = 0; ihit < nBgoHits; ++ihit)
	{
		// Get layer ID
		auto layerID = bgohits->GetLayerID(ihit);
		// Get bar global ID
		auto iBar = ((bgohits->fGlobalBarID)[ihit] >> 6) & 0x1f;
		layerBarIndex[layerID].push_back(ihit);
		layerBarNumber[layerID].push_back(iBar);
	}

	for (int lay = 0; lay < nLayers; ++lay)
	{
		// Setting default value for maximum bar index and energy for each layer
		int max_hit 			{-1};			// hit of the maximum energy release in a layer
		size_t max_hit_idx 		{0};			// hit index relative to the maximum energy release in a layer	
		double maxE 			{0};			// maximum energy release in a layer
	
		// Find the maximum of the nergy release in a certain layer, together with the bar ID
		for (size_t layerBar_idx=0; layerBar_idx<layerBarIndex[lay].size(); ++layerBar_idx)
		{	
			double hitE = (bgohits->fEnergy)[layerBarIndex[lay][layerBar_idx]];
			layerBarEnergy[lay][layerBarNumber[lay][layerBar_idx]] = hitE;
			if (hitE > maxE)
			{
				maxE = hitE;
				max_hit = layerBarIndex[lay][layerBar_idx];
				max_hit_idx = layerBar_idx;
			}
		}

		iMaxLayer[lay] 	= max_hit;
		rmsLayer[lay] 	= 0;
		fracLayer[lay] 	= 0;
		eLayer[lay] 	= 0;

		if (maxE > layerMinEnergy)
		{
			// Find the bar index regarding the maximum energy release in a certain layer
			idxBarMaxLayer[lay] = layerBarNumber[lay][max_hit_idx];
			// Register the maximum energy release of a layer
			eCoreLayer[lay] = maxE;
			// Find the coordinate (weighted by the energy release) of the bar with the biggest energy release in a certain layer
			double coordMax = lay % 2 ? bgohits->GetHitX(max_hit) : bgohits->GetHitY(max_hit);
			eCoreCoord[lay] = maxE * coordMax;
			// Consider the nearest bar respect to the max one in order to better interpolate the position
			if (max_hit_idx > 0 && max_hit_idx < 21)
			{
				for (auto it = std::begin(layerBarIndex[lay]); it != std::end(layerBarIndex[lay]); ++it)
				{
					auto iBar = ((bgohits->fGlobalBarID)[*it] >> 6) & 0x1f;
					if (fabs(iBar - max_hit_idx) == 1)
					{
						double hitE = (bgohits->fEnergy)[*it];
						double thisCoord = lay % 2 ? bgohits->GetHitX(*it) : bgohits->GetHitY(*it);
						eCoreLayer[lay] += hitE;
						eCoreCoord[lay] += hitE * thisCoord;
					}
				}
			}
			// Get the CoG coordinate of the max energy bar cluster
			eCoreCoord[lay] /= eCoreLayer[lay];
			// Get the energy RMS of a layer
			for (auto it = std::begin(layerBarIndex[lay]); it != std::end(layerBarIndex[lay]); ++it)
			{
				auto hitE = (bgohits->fEnergy)[*it];
				auto thisCoord = lay % 2 ? bgohits->GetHitX(*it) : bgohits->GetHitY(*it);
				eLayer[lay] += hitE;
				rmsLayer[lay] += hitE * pow((thisCoord - eCoreCoord[lay]), 2);
			}
			rmsLayer[lay] = sqrt(rmsLayer[lay] / eLayer[lay]);
			fracLayer[lay] = eLayer[lay] / bgoTotalE;
			if (layerBarNumber[lay].size() <= 1)
				rmsLayer[lay] = 0;
		}
		sumRms += rmsLayer[lay];
	}

	// Find last energy layer
	if (lastLayer == -999)
		for (auto it = std::rbegin(fracLayer); it != std::rend(fracLayer); ++it)
			if (*it)
			{
				lastLayer = fracLayer.size() - 1 - std::distance(fracLayer.rbegin(), it);
				break;
			}

	// Calculate rvalue
	calculate_rvalue();
	// Calculate lvalue
	calculate_lvalue(bgoTotalE);
}

double shower_profile_fit_function(double *x, double *par)
{
    double xx = x[0];
    double func = par[0] * par[1] * (TMath::Power(par[1] * xx, par[2] - 1) * TMath::Exp(-par[1] * xx)) / TMath::Gamma(par[2]);
    return func;
}

void DmpBgoContainer::calculate_lvalue(const double bgoTotalE)
{
	const double bgoX0 				{11.2};     // BGO X0 in mm
    const double bgoEc 				{10.50};    // BGO critical energy for electrons in MeV
    const double b_shower_par 		{0.5}; 		// b parameter electromagnetic shower in BGO
   
	double a_shower_par {1 + b_shower_par * (TMath::Log(bgoTotalE/bgoEc) - 0.5)}; 	// a parameter electromagnetic shower in BGO

	for (int idx = 0; idx < DAMPE_bgo_nLayers; ++idx)
        t_bgo[idx] = t_bgo_norm[idx] = (BGO_bar_lateral * (idx + 1)) / (bgoX0 * trajectoryDirection2D.CosTheta());
	
	// normalize t_bgo respect to the maximum value. This way values should be in the range [0,1]
	double max_t_bgo_no_cosine_correction {(BGO_bar_lateral * (DAMPE_bgo_nLayers + 1))/bgoX0};
	for (int idx = 0; idx < DAMPE_bgo_nLayers; ++idx)
		//t_bgo_norm[idx] /= t_bgo_norm[DAMPE_bgo_nLayers-1];
		t_bgo_norm[idx] /= max_t_bgo_no_cosine_correction;

	// Build the TGraph
	TGraph shower_profile_gr(DAMPE_bgo_nLayers, &t_bgo[0], &eLayer[0]);
	TGraph shower_profile_gr_norm(DAMPE_bgo_nLayers, &t_bgo_norm[0], &eLayer[0]);
	// Build the fit function
	TF1 fitfunc("fitfunc", shower_profile_fit_function, t_bgo[0], t_bgo[DAMPE_bgo_nLayers-1], 3);
	TF1 fitfunc_norm("fitfunc_norm", shower_profile_fit_function, t_bgo_norm[0], t_bgo_norm[DAMPE_bgo_nLayers-1], 3);
	fitfunc.SetParameter(0, bgoTotalE);
	fitfunc.SetParameter(1, b_shower_par);
	fitfunc.SetParameter(2, a_shower_par);
	fitfunc_norm.SetParameter(0, bgoTotalE);
	fitfunc_norm.SetParameter(1, b_shower_par);
	fitfunc_norm.SetParameter(2, a_shower_par);
	// Fit the shower profile
	shower_profile_gr.Fit(&fitfunc, "qR");
	shower_profile_gr_norm.Fit(&fitfunc_norm, "qR");
	// Extract the X value corresponding to the maximum value of the fit function
	maximum_shower_position = fitfunc.GetMaximumX();
	maximum_shower_position_norm = fitfunc_norm.GetMaximumX();
	// Get the lvalue
	//lvalue = maximum_shower_position_norm;
	//lvalue = (maximum_shower_position - t_bgo[0])/(t_bgo[DAMPE_bgo_nLayers-1] - t_bgo[0]);
	lvalue = (maximum_shower_position - t_bgo[0])/(max_t_bgo_no_cosine_correction - t_bgo[0]);
}

void DmpBgoContainer::calculate_rvalue()
{
	rvalue = 0;
	for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
	{
		if (iMaxLayer[lIdx] > -1)
		{
			if (idxBarMaxLayer[lIdx] != 0 && idxBarMaxLayer[lIdx] != 21)
			{
				if (eLayer[lIdx])
					rvalue += layerBarEnergy[lIdx][idxBarMaxLayer[lIdx]]/eLayer[lIdx];
			}
			else
				break;
		}
	}
}

std::vector<short> DmpBgoContainer::GetSingleLayerBarNumber(int nLayer)
{
	return layerBarNumber[nLayer];
}

std::vector<short> DmpBgoContainer::GetSingleLayerBarIndex(int nLayer) 
{
	return layerBarIndex[nLayer];
}

std::vector<std::vector<short>> DmpBgoContainer::GetLayerBarNumber()
{
	return layerBarNumber;
}

std::vector<double> DmpBgoContainer::GetRmsLayer()
{
	return rmsLayer;
}

const double DmpBgoContainer::GetSumRMS()
{
	return sumRms;
}

std::vector<double> DmpBgoContainer::GetFracLayer()
{
	return fracLayer;
}

const double DmpBgoContainer::GetSingleFracLayer(const int lIdx)
{

	if (lIdx >= 0 && (unsigned int)lIdx <= fracLayer.size() - 1)
		return fracLayer[lIdx];
	else
		return -1;
}

const int DmpBgoContainer::GetLastEnergyLayer()
{
	return lastLayer;
}

std::vector<int> DmpBgoContainer::GetIdxBarMaxLayer()
{
	return idxBarMaxLayer;
}

const int DmpBgoContainer::GetIdxBarMaxSingleLayer(const int layedIdx)
{
	return idxBarMaxLayer[layedIdx];
}

std::vector<double> DmpBgoContainer::GetELayer()
{
	return eLayer;
}

const double DmpBgoContainer::GetESingleLayer(const int layedIdx)
{
	return eLayer[layedIdx];
}

std::vector<int> DmpBgoContainer::GetiMaxLayer()
{
	return iMaxLayer;
}

const int DmpBgoContainer::GetiMaxSingleLayer(const int layedIdx)
{
	return iMaxLayer[layedIdx];
}

const int DmpBgoContainer::GetNhits()
{
	return nBgoHits;
}

const std::vector<double> DmpBgoContainer::GetLayerEnergies()
{
	return eLayer;
}

const std::vector<std::vector<double>> DmpBgoContainer::GetLayerBarEnergies()
{
	return layerBarEnergy;
}

const std::vector<double> DmpBgoContainer::GetBGOslope()
{
	return slope;
}

const std::vector<double> DmpBgoContainer::GetBGOintercept()
{
	return intercept;
}

const std::vector<double> DmpBgoContainer::FastBGOslope(const std::shared_ptr<DmpEvtBgoRec> bgorec)
{
	slope[0] = bgorec->GetSlopeXZ();
	slope[1] = bgorec->GetSlopeYZ();
	return slope;
}

const std::vector<double> DmpBgoContainer::FastBGOintercept(const std::shared_ptr<DmpEvtBgoRec> bgorec)
{
	intercept[0] = bgorec->GetInterceptXZ();
	intercept[1] = bgorec->GetInterceptYZ();
	return intercept;
}

const TVector3 DmpBgoContainer::GetBGOTrajectory2D()
{
	return trajectoryDirection2D;
}

const double DmpBgoContainer::GetRValue()
{
	return rvalue;
}

const double DmpBgoContainer::GetLValue()
{
	return lvalue;
}

const double DmpBgoContainer::GetMaximumShowerPosition()
{
	return maximum_shower_position;
}

const double DmpBgoContainer::GetMaximumShowerPositionNorm()
{
	return maximum_shower_position_norm;
}

const std::vector<double> DmpBgoContainer::GetTShowerProfile()
{
	return t_bgo;
}

const std::vector<double> DmpBgoContainer::GetTShowerProfileNorm()
{
	return t_bgo_norm;
}

void DmpBgoContainer::Reset()
{
	layerBarIndex 					= std::vector<std::vector<short>>(DAMPE_bgo_nLayers, std::vector<short>());
	layerBarNumber 					= std::vector<std::vector<short>>(DAMPE_bgo_nLayers, std::vector<short>());
	layerBarEnergy 					= std::vector<std::vector<double>>(DAMPE_bgo_nLayers, std::vector<double> (DAMPE_bgo_bars_layer, -999));
	idxBarMaxLayer 					= std::vector<int>(DAMPE_bgo_nLayers, -1);
	iMaxLayer						= std::vector<int>(DAMPE_bgo_nLayers, -1);
	rmsLayer 						= std::vector<double>(DAMPE_bgo_nLayers, 0);
	fracLayer						= std::vector<double>(DAMPE_bgo_nLayers, 0);
	eLayer							= std::vector<double>(DAMPE_bgo_nLayers, 0);
	eCoreLayer						= std::vector<double>(DAMPE_bgo_nLayers, 0);
	eCoreCoord						= std::vector<double>(DAMPE_bgo_nLayers, 0);
	slope							= std::vector<double>(2, -999);
	intercept						= std::vector<double>(2, -999);

	nBgoHits						= -999;
	lastLayer						= -999;
	sumRms							= 0;

	trajectoryDirection2D 			= TVector3();

	rvalue							= 0;
	lvalue							= -999;
	maximum_shower_position 		= -999;
	maximum_shower_position_norm 	= -999;
	t_bgo							= std::vector<double>(DAMPE_bgo_nLayers, 0);
	t_bgo_norm						= std::vector<double>(DAMPE_bgo_nLayers, 0);
}