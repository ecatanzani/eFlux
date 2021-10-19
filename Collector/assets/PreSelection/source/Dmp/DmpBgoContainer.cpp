#include "Dmp/DmpBgoContainer.h"

#include <cmath>
#include <iostream>

void DmpBgoContainer::scanBGOHits(
	const std::shared_ptr<DmpEvtBgoHits> bgohits,
	const std::shared_ptr<DmpEvtBgoRec> bgorec,
	const double bgoTotalE,
	const double layerMinEnergy,
	const int nLayers)
{
	// Compute BGO shower axis
	if (slope[0] == -999 && slope[1] == -999) {
		slope[0] = bgorec->GetSlopeXZ();
		slope[1] = bgorec->GetSlopeYZ();
	}
	if (intercept[0] == -999 && intercept[1] == -999) {
		intercept[0] = bgorec->GetInterceptXZ();
		intercept[1] = bgorec->GetInterceptYZ();
	}
	// Get BGO reco 2D trajectory
	trajectoryDirection2D = bgorec->GetTrajectoryDirection2D();
	
	// Get the number of BGO hits
	nBgoHits = bgohits->GetHittedBarNumber();

	// Vector of BGO Z values
	std::vector<double> bgolayersZ(DAMPE_bgo_nLayers, -999);
	// Scan BGO hits
	for (int ihit = 0; ihit < nBgoHits; ++ihit) {
		// Get layer ID
		auto layerID = bgohits->GetLayerID(ihit);
		if (bgolayersZ[layerID] == -999)
			bgolayersZ[layerID] = BGO_TopZ + BGO_layer_displacement * layerID;
		// Get bar global ID
		auto iBar = ((bgohits->fGlobalBarID)[ihit] >> 6) & 0x1f;
		layerBarIndex[layerID].push_back(ihit);
		layerBarNumber[layerID].push_back(iBar);
	}

	for (int lay = 0; lay < nLayers; ++lay) {
		// Setting default value for maximum bar index and energy for each layer
		int imax {-1};
		double maxE {0};
	
		// Find the maximum of the nergy release in a certain layer, together with the bar ID
		for (auto it = std::begin(layerBarNumber[lay]); it != std::end(layerBarNumber[lay]); ++it) {
			auto index = std::distance(std::begin(layerBarNumber[lay]), it);
			int ihit = layerBarIndex[lay][index];
			double hitE = (bgohits->fEnergy)[ihit];
			layerBarEnergy[lay][index] = hitE;
			if (hitE > maxE) {
				maxE = hitE;
				imax = ihit;
			}
		}

		iMaxLayer[lay] = imax;
		rmsLayer[lay] = 0;
		fracLayer[lay] = 0;
		eLayer[lay] = 0;

		if (maxE > layerMinEnergy) {
			// Find the bar index regarding the maximum energy release in a certain layer
			auto iBarMax = ((bgohits->fGlobalBarID)[imax] >> 6) & 0x1f;
			idxBarMaxLayer[lay] = iBarMax;
			// Register the maximum energy release of a layer
			eCoreLayer[lay] = maxE;
			// Find the coordinate (weighted by the nergy release) of the bar with the biggest energy release in a certain layer
			double coordMax {lay % 2 ? bgohits->GetHitX(imax) : bgohits->GetHitY(imax)};
			eCoreCoord[lay] = maxE * coordMax;
			// Consider the nearest bar respect to the max one in order to better interpolate the position
			if (iBarMax > 0 && iBarMax < 21) {
				for (auto it = std::begin(layerBarNumber[lay]); it != std::end(layerBarNumber[lay]); ++it) {
					auto index = std::distance(std::begin(layerBarNumber[lay]), it);
					int ihit = layerBarIndex[lay][index];
					auto iBar = ((bgohits->fGlobalBarID)[ihit] >> 6) & 0x1f;
					if (abs(iBar - iBarMax) == 1) {
					//if (iBar - iBarMax == 1 || iBar - iBarMax == -1) {
						double hitE = (bgohits->fEnergy)[ihit];
						double thisCoord {lay % 2 ? bgohits->GetHitX(ihit) : bgohits->GetHitY(ihit)};
						eCoreLayer[lay] += hitE;
						eCoreCoord[lay] += hitE * thisCoord;
					}
				}
			}
			// Get the CoG coordinate of the max energy bar cluster
			eCoreCoord[lay] /= eCoreLayer[lay];
			// Get the energy RMS of a layer
			for (auto it = std::begin(layerBarNumber[lay]); it != std::end(layerBarNumber[lay]); ++it) {
				auto index = std::distance(std::begin(layerBarNumber[lay]), it);
				int ihit = layerBarIndex[lay][index];
				auto hitE = (bgohits->fEnergy)[ihit];
				auto thisCoord {lay % 2 ? bgohits->GetHitX(ihit) : bgohits->GetHitY(ihit)};
				eLayer[lay] += hitE;
				rmsLayer[lay] += hitE * pow(thisCoord - eCoreCoord[lay], 2);
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
			if (*it) {
				lastLayer = fracLayer.size() - 1 - std::distance(std::rbegin(fracLayer), it);
				break;
			}
}

std::vector<short> DmpBgoContainer::GetSingleLayerBarNumber(int nLayer) {
	return layerBarNumber[nLayer];
}

std::vector<std::vector<short>> DmpBgoContainer::GetLayerBarNumber() {
	return layerBarNumber;
}

std::vector<double> DmpBgoContainer::GetRmsLayer() {
	return rmsLayer;
}

const double DmpBgoContainer::GetSumRMS() {
	return sumRms;
}

std::vector<double> DmpBgoContainer::GetFracLayer() {
	return fracLayer;
}

const double DmpBgoContainer::GetSingleFracLayer(const int lIdx) {

	if (lIdx >= 0 && static_cast<unsigned int>(lIdx) <= fracLayer.size() - 1)
		return fracLayer[lIdx];
	else
		return -1;
}

const int DmpBgoContainer::GetLastEnergyLayer() {
	return lastLayer;
}

std::vector<int> DmpBgoContainer::GetIdxBarMaxLayer() {
	return idxBarMaxLayer;
}

const int DmpBgoContainer::GetIdxBarMaxSingleLayer(const int layedIdx) {
	return idxBarMaxLayer[layedIdx];
}

std::vector<double> DmpBgoContainer::GetELayer() {
	return eLayer;
}

const double DmpBgoContainer::GetESingleLayer(const int layedIdx) {
	return eLayer[layedIdx];
}

std::vector<int> DmpBgoContainer::GetiMaxLayer() {
	return iMaxLayer;
}

const int DmpBgoContainer::GetiMaxSingleLayer(const int layedIdx) {
	return iMaxLayer[layedIdx];
}

const int DmpBgoContainer::GetNhits() {
	return nBgoHits;
}

const std::vector<double> DmpBgoContainer::GetLayerEnergies() {
	return eLayer;
}

const std::vector<std::vector<double>> DmpBgoContainer::GetLayerBarEnergies() {
	return layerBarEnergy;
}

const std::vector<double> DmpBgoContainer::GetBGOslope() {
	return slope;
}

const std::vector<double> DmpBgoContainer::GetBGOintercept() {
	return intercept;
}

const std::vector<double> DmpBgoContainer::FastBGOslope(const std::shared_ptr<DmpEvtBgoRec> bgorec) {
	slope[0] = bgorec->GetSlopeXZ();
	slope[1] = bgorec->GetSlopeYZ();
	return slope;
}

const std::vector<double> DmpBgoContainer::FastBGOintercept(const std::shared_ptr<DmpEvtBgoRec> bgorec) {
	intercept[0] = bgorec->GetInterceptXZ();
	intercept[1] = bgorec->GetInterceptYZ();
	return intercept;
}

const TVector3 DmpBgoContainer::GetBGOTrajectory2D() {
	return trajectoryDirection2D;
}