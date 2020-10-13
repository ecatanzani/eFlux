#ifndef DMPBGOCONTAINER_H
#define DMPBGOCONTAINER_H

#include <vector>
#include <memory>

#include "data_cuts.h"
#include "DmpEvtBgoHits.h"
#include "DAMPE_geo_structure.h"

class DmpBgoContainer
{
public:
	DmpBgoContainer()	:	layerBarIndex(DAMPE_bgo_nLayers, std::vector<short>()),
							layerBarNumber(DAMPE_bgo_nLayers, std::vector<short>()),
							idxBarMaxLayer(DAMPE_bgo_nLayers, -1),
							iMaxLayer(DAMPE_bgo_nLayers, -1),
							rmsLayer(DAMPE_bgo_nLayers, 0),
							fracLayer(DAMPE_bgo_nLayers, 0),
							eLayer(DAMPE_bgo_nLayers, 0),
							eCoreLayer(DAMPE_bgo_nLayers, 0),
							eCoreCoord(DAMPE_bgo_nLayers, 0)
							{
							}

	DmpBgoContainer(int m_size)		:	layerBarIndex(m_size, std::vector<short>()),
								  		layerBarNumber(m_size, std::vector<short>()),
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
		const double bgoTotalE,
		const cuts_conf &data_cuts,
		const int nLayers = DAMPE_bgo_nLayers);
	std::vector<short> GetSingleLayerBarNumber(int nLayer);
	std::vector<std::vector<short>> GetLayerBarNumber();
	std::vector<double> GetRmsLayer();
	std::vector<double> GetFracLayer();
	const double GetSingleFracLayer(const int lIdx);
	const double GetLastFracLayer();
	const double GetLastFFracLayer();
	std::vector<int> GetIdxBarMaxLayer();
	const int GetIdxLastLayer();
	const int GetFracIdxLastLayer();
	const int GetIdxBarMaxSingleLayer(const int layedIdx);
	std::vector<double> GetELayer();
	const double GetESingleLayer(const int layedIdx);
	std::vector<int> GetiMaxLayer();
	const int GetiMaxSingleLayer(const int layedIdx);
	const double GetSumRMS();
	const int GetNhits();

private:
	std::vector<std::vector<short>> layerBarIndex;	// arrange BGO hits by layer
	std::vector<std::vector<short>> layerBarNumber; // arrange BGO bars by layer
	std::vector<int> idxBarMaxLayer;
	std::vector<double> rmsLayer;
	std::vector<double> fracLayer;
	std::vector<double> eLayer;
	std::vector<double> eCoreLayer;
	std::vector<double> eCoreCoord;
	std::vector<int> iMaxLayer;
	int nBgoHits = 0;
	int lastLayer = 0;
	double sumRms = 0;
};

#endif