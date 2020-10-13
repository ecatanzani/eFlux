#ifndef DMPPSDCONTAINER_H
#define DMPPSDCONTAINER_H

#include <vector>
#include <memory>

#include "data_cuts.h"
#include "DmpEvtPsdHits.h"
#include "DAMPE_geo_structure.h"

class DmpPsdContainer
{
    public:
        DmpPsdContainer()   :   layerBarIndexPsd(DAMPE_bgo_nLayers, std::vector<short>()),
                                layerBarNumberPsd(DAMPE_bgo_nLayers, std::vector<short>()),
                                layerBarEnergyPsd(DAMPE_bgo_nLayers, std::vector<double>()),
                                layerBarUsedPsd(DAMPE_bgo_nLayers, std::vector<short>()),
                                psdCluster_idxBeg(DAMPE_bgo_nLayers, std::vector<short>()),
                                psdCluster_length(DAMPE_bgo_nLayers, std::vector<short>()),
                                psdCluster_idxMaxE(DAMPE_bgo_nLayers, std::vector<short>()),
                                psdCluster_E(DAMPE_bgo_nLayers, std::vector<double>()),
                                psdCluster_maxE(DAMPE_bgo_nLayers, std::vector<double>()),
                                psdCluster_maxEcoordinate(DAMPE_bgo_nLayers, std::vector<double>()),
                                psdCluster_coordinate(DAMPE_bgo_nLayers, std::vector<double>()),
                                psdCluster_Z(DAMPE_bgo_nLayers, std::vector<double>())
                                {
                                }

        DmpPsdContainer(int msize)  :   layerBarIndexPsd(msize, std::vector<short>()),
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
            const cuts_conf &data_cuts,
            const int nLayers = DAMPE_psd_nLayers);
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
        unsigned int nPsdClusters;
        unsigned int nPsdClustersX;
        unsigned int nPsdClustersY;
};

#endif