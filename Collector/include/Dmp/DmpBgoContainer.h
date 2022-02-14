#ifndef DMPBGOCONTAINER_H
#define DMPBGOCONTAINER_H

#include <vector>
#include <memory>

#include "DmpEvtBgoRec.h"
#include "DmpEvtBgoHits.h"
#include "Dmp/DmpGeoStruct.h"

#include "TVector3.h"

class DmpBgoContainer
{
public:
	DmpBgoContainer() : layerBarIndex(DAMPE_bgo_nLayers, std::vector<short>()),
						layerBarNumber(DAMPE_bgo_nLayers, std::vector<short>()),
						layerBarEnergy(DAMPE_bgo_nLayers, std::vector<double> (DAMPE_bgo_bars_layer, -999)),
						idxBarMaxLayer(DAMPE_bgo_nLayers, -1),
						iMaxLayer(DAMPE_bgo_nLayers, -1),
						rmsLayer(DAMPE_bgo_nLayers, 0),
						fracLayer(DAMPE_bgo_nLayers, 0),
						eLayer(DAMPE_bgo_nLayers, 0),
						eCoreLayer(DAMPE_bgo_nLayers, 0),
						eCoreCoord(DAMPE_bgo_nLayers, 0)
	{
	}

	DmpBgoContainer(int m_size) : layerBarIndex(m_size, std::vector<short>()),
								  layerBarNumber(m_size, std::vector<short>()),
								  layerBarEnergy(m_size, std::vector<double> (DAMPE_bgo_bars_layer, -999)),
								  idxBarMaxLayer(m_size, -1),
								  iMaxLayer(m_size, -1),
								  rmsLayer(m_size, 0),
								  fracLayer(m_size, 0),
								  eLayer(m_size, 0),
								  eCoreLayer(m_size, 0),
								  eCoreCoord(m_size, 0)
	{
	}

	~DmpBgoContainer(){};
	void scanBGOHits(
		const std::shared_ptr<DmpEvtBgoHits> bgohits,
		const std::shared_ptr<DmpEvtBgoRec> bgorec,
		const double bgoTotalE,
		const double layerMinEnergy,
		const int nLayers = DAMPE_bgo_nLayers);

	// Reset
	void Reset();
	
	// BGO Bars
	std::vector<short> GetSingleLayerBarNumber(int nLayer);
	std::vector<short> GetSingleLayerBarIndex(int nLayer);
	std::vector<std::vector<short>> GetLayerBarNumber();
	std::vector<int> GetIdxBarMaxLayer();
	const int GetIdxBarMaxSingleLayer(const int layedIdx);
	// BGO RMS
	std::vector<double> GetRmsLayer();
	const double GetSumRMS();
	// BGO energy fraction
	std::vector<double> GetFracLayer();
	const double GetSingleFracLayer(const int lIdx);
	// BGO last energy layer
	const int GetLastEnergyLayer();
	// BGO energy
	std::vector<double> GetELayer();
	const double GetESingleLayer(const int layedIdx);
	std::vector<int> GetiMaxLayer();
	const int GetiMaxSingleLayer(const int layedIdx);
	const std::vector<double> GetLayerEnergies();
	const std::vector<std::vector<double>> GetLayerBarEnergies();
	// BGO slope and itercept
	const std::vector<double> GetBGOslope();
	const std::vector<double> GetBGOintercept();
	const std::vector<double> FastBGOslope(const std::shared_ptr<DmpEvtBgoRec> bgorec);
	const std::vector<double> FastBGOintercept(const std::shared_ptr<DmpEvtBgoRec> bgorec);
	const TVector3 GetBGOTrajectory2D();
	// BGO hits
	const int GetNhits();

private:
	// BGO variables
	std::vector<std::vector<short>> layerBarIndex;	// arrange BGO hits by layer
	std::vector<std::vector<short>> layerBarNumber; // arrange BGO bars by layer
	std::vector<std::vector<double>> layerBarEnergy; // store BGO bars energy
	std::vector<int> idxBarMaxLayer;
	std::vector<double> rmsLayer;
	std::vector<double> fracLayer;
	std::vector<double> eLayer;
	std::vector<double> eCoreLayer;
	std::vector<double> eCoreCoord;
	std::vector<int> iMaxLayer;
	int nBgoHits = -999;
	int lastLayer = -999;
	double sumRms = 0;
	std::vector<double> slope{-999, -999};
	std::vector<double> intercept{-999, -999};
	TVector3 trajectoryDirection2D;
};

#endif