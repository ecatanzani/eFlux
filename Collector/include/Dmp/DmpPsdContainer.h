#ifndef DMPPSDCONTAINER_H
#define DMPPSDCONTAINER_H

#include <vector>
#include <memory>

#include "DmpEvtPsdHits.h"
#include "Dmp/DmpGeoStruct.h"

class DmpPsdContainer
{
public:
	DmpPsdContainer() : layerBarIndexPsd(DAMPE_psd_nLayers, std::vector<short>()),
						layerBarNumberPsd(DAMPE_psd_nLayers, std::vector<short>()),
						layerBarEnergyPsd(DAMPE_psd_nLayers, std::vector<double>()),
						layerBarUsedPsd(DAMPE_psd_nLayers, std::vector<short>()),
						psdCluster_idxBeg(DAMPE_psd_nLayers, std::vector<short>()),
						psdCluster_length(DAMPE_psd_nLayers, std::vector<short>()),
						psdCluster_idxMaxE(DAMPE_psd_nLayers, std::vector<short>()),
						psdCluster_E(DAMPE_psd_nLayers, std::vector<double>()),
						psdCluster_maxE(DAMPE_psd_nLayers, std::vector<double>()),
						psdCluster_maxEcoordinate(DAMPE_psd_nLayers, std::vector<double>()),
						psdCluster_coordinate(DAMPE_psd_nLayers, std::vector<double>()),
						psdCluster_Z(DAMPE_psd_nLayers, std::vector<double>())
	{
	}

	DmpPsdContainer(int msize) : layerBarIndexPsd(msize, std::vector<short>()),
								 layerBarNumberPsd(msize, std::vector<short>()),
								 layerBarEnergyPsd(msize, std::vector<double>()),
								 layerBarUsedPsd(msize, std::vector<short>()),
								 psdCluster_idxBeg(msize, std::vector<short>()),
								 psdCluster_length(msize, std::vector<short>()),
								 psdCluster_idxMaxE(msize, std::vector<short>()),
								 psdCluster_E(msize, std::vector<double>()),
								 psdCluster_maxE(msize, std::vector<double>()),
								 psdCluster_maxEcoordinate(msize, std::vector<double>()),
								 psdCluster_coordinate(msize, std::vector<double>()),
								 psdCluster_Z(msize, std::vector<double>())
	{
	}

	~DmpPsdContainer(){};
	void scanPSDHits(
		const std::shared_ptr<DmpEvtPsdHits> psdhits,
		const double PSD_bar_min_energy_release,
		const int nLayers = DAMPE_psd_nLayers);
	void Reset();
	std::vector<std::vector<short>> getPsdClusterIdxBegin();
	std::vector<std::vector<double>> getPsdClusterZ();
	std::vector<std::vector<double>> getPsdClusterMaxE();
	std::vector<std::vector<short>> getPsdClusterIdxMaxE();
	std::vector<std::vector<double>> getPsdClusterMaxECoo();
	std::vector<short> getGlobalBarID();
	std::vector<double> getHitZ();
	const unsigned int getPsdNclusters();
	const unsigned int getPsdNclustersX();
	const unsigned int getPsdNclustersY();

private:
	std::vector<std::vector<short>> layerBarIndexPsd;	// Arrange PSD hit index
	std::vector<std::vector<short>> layerBarNumberPsd;	// Arrange PSD unique bar number
	std::vector<std::vector<double>> layerBarEnergyPsd; // Arrange PSD hit energy
	std::vector<std::vector<short>> layerBarUsedPsd;	// Mark PSD hits index used for clustering

	// PSD clusters
	std::vector<std::vector<short>> psdCluster_idxBeg;
	std::vector<std::vector<short>> psdCluster_length;
	std::vector<std::vector<short>> psdCluster_idxMaxE;
	std::vector<std::vector<double>> psdCluster_E;
	std::vector<std::vector<double>> psdCluster_maxE;
	std::vector<std::vector<double>> psdCluster_maxEcoordinate; // Arrange X and Y coordinates regarding the max energy release on PSD layer
	std::vector<std::vector<double>> psdCluster_coordinate;		// weighted if more than 1 strip
	std::vector<std::vector<double>> psdCluster_Z;				// weighted if more than 1 strip

	// Z position of PSD hits
	std::vector<double> hitZ;

	// Global bar ID
	std::vector<short> globalBarID;

	// PSD clusters
	unsigned int nPsdClusters 		{0};	
	unsigned int nPsdClustersX		{0};
	unsigned int nPsdClustersY		{0};
};

#endif