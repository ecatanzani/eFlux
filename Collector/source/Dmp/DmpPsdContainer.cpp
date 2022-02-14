#include "Dmp/DmpPsdContainer.h"

void DmpPsdContainer::scanPSDHits(
	const std::shared_ptr<DmpEvtPsdHits> psdhits,
	const double PSD_bar_min_energy_release,
	const int nLayers)
{
	// Find max index and energy release for PSD hits
	std::vector<int> iMaxBarPsd(2, -1);
	std::vector<double> eMaxBarPsd(2, 0);
	bool first_hit = true;

	// Loop on PSD hits
	int nPSD_tot_entries = psdhits->GetHittedBarNumber();
	hitZ.resize(nPSD_tot_entries);
	globalBarID = static_cast<std::vector<short>>(psdhits->fGlobalBarID);
	for (int ihit = 0; ihit < nPSD_tot_entries; ++ihit)
	{
		double hitE = (psdhits->fEnergy)[ihit];
		hitZ[ihit] = psdhits->GetHitZ(ihit);
		if (hitE > 0.5)
		{
			short layerID = psdhits->GetLayerID(ihit);
			int iBar = ((psdhits->fGlobalBarID)[ihit] >> 6);
			layerBarUsedPsd[layerID].push_back(0);
			layerBarIndexPsd[layerID].push_back(ihit);
			layerBarNumberPsd[layerID].push_back(iBar);
			layerBarEnergyPsd[layerID].push_back(hitE);
			if (first_hit)
			{
				eMaxBarPsd[layerID] = hitE;
				iMaxBarPsd[layerID] = layerBarIndexPsd[layerID].size() - 1;
				first_hit = false;
			}
			else
			{
				if (hitE > eMaxBarPsd[layerID])
				{
					eMaxBarPsd[layerID] = hitE;
					iMaxBarPsd[layerID] = layerBarIndexPsd[layerID].size() - 1;
				}
			}
		}
	}

	/*
		cluster finding : starting from maxE, seed + the highest E neigboring bar
		coordinate: weighted average of two highest neigboring bars
	*/

	for (int nLayer = 0; nLayer < nLayers; ++nLayer)
	{
		while (eMaxBarPsd[nLayer] > PSD_bar_min_energy_release)
		{
			int psdhits_maxIdx = layerBarIndexPsd[nLayer][iMaxBarPsd[nLayer]];
			int psdbar_maxIdx = layerBarNumberPsd[nLayer][iMaxBarPsd[nLayer]];
			double cluster_energy = eMaxBarPsd[nLayer];
			double actual_coordinate = nLayer % 2 ? psdhits->GetHitX(psdhits_maxIdx) : psdhits->GetHitY(psdhits_maxIdx);
			psdCluster_maxEcoordinate[nLayer].push_back(actual_coordinate);
			double actualZ = psdhits->GetHitZ(psdhits_maxIdx);
			layerBarUsedPsd[nLayer][iMaxBarPsd[nLayer]] = 1;

			int cluster_closest_barIdx = -1;
			double energy_leftMaxBar = 0;
			double energy_rightMaxBar = 0;
			int cluster_fIdx = psdhits_maxIdx; // First index of the cluster
			int cluster_size = 1;			   // Size of the cluster

			// Complete the cluster building procedure with the closest bars to the one with the maximum energy release
			if ((iMaxBarPsd[nLayer] - 1) > 0) // Check the previous bar
				if ((layerBarUsedPsd[nLayer][iMaxBarPsd[nLayer] - 1] == 0) && (psdbar_maxIdx == layerBarNumberPsd[nLayer][iMaxBarPsd[nLayer] - 1] + 1))
					energy_leftMaxBar = layerBarEnergyPsd[nLayer][iMaxBarPsd[nLayer] - 1];
			if ((unsigned int)(iMaxBarPsd[nLayer] + 1) < layerBarEnergyPsd[nLayer].size()) // Check the next bar
				if ((layerBarUsedPsd[nLayer][iMaxBarPsd[nLayer] + 1] == 0) && (psdbar_maxIdx == layerBarNumberPsd[nLayer][iMaxBarPsd[nLayer] + 1] - 1))
					energy_rightMaxBar = layerBarEnergyPsd[nLayer][iMaxBarPsd[nLayer] + 1];

			/* 
				FINAL CLUSTER BUILING:
				Now we have the central bar with the highest energy release, the left and the right one. 
				One cluster is composed by 2 bars. We have to choose, between the 2 lateral bars, that one with the highest energy release
			*/

			if (energy_leftMaxBar > energy_rightMaxBar)
			{
				cluster_closest_barIdx = psdhits_maxIdx - 1;
				cluster_energy += energy_leftMaxBar;
				++cluster_size;
				cluster_fIdx = psdhits_maxIdx - 1;
				layerBarUsedPsd[nLayer][iMaxBarPsd[nLayer] - 1] = 1;
			}
			//else if (energy_leftMaxBar < energy_rightMaxBar)
			else
			{
				cluster_closest_barIdx = psdhits_maxIdx + 1;
				cluster_energy += energy_rightMaxBar;
				++cluster_size;
				layerBarUsedPsd[nLayer][iMaxBarPsd[nLayer] + 1] = 1;
			}

			if (cluster_closest_barIdx > -1)
			{
				double cluster_closest_bar_coordinate = iMaxBarPsd[nLayer] % 2 ? psdhits->GetHitX(cluster_closest_barIdx) : psdhits->GetHitY(cluster_closest_barIdx);
				double cluster_closest_energy = (psdhits->fEnergy)[cluster_closest_barIdx];
				double cluster_closestZ = psdhits->GetHitZ(cluster_closest_barIdx);
				actual_coordinate = (actual_coordinate * eMaxBarPsd[nLayer] + cluster_closest_bar_coordinate * cluster_closest_energy) / (eMaxBarPsd[nLayer] + cluster_closest_energy);
				actualZ = (actualZ * eMaxBarPsd[nLayer] + cluster_closestZ * cluster_closest_energy) / (eMaxBarPsd[nLayer] + cluster_closest_energy);
			}

			psdCluster_idxBeg[nLayer].push_back(cluster_fIdx);
			psdCluster_length[nLayer].push_back(cluster_size);
			psdCluster_idxMaxE[nLayer].push_back(psdhits_maxIdx);
			psdCluster_maxE[nLayer].push_back(eMaxBarPsd[nLayer]);
			psdCluster_E[nLayer].push_back(cluster_energy);
			psdCluster_coordinate[nLayer].push_back(actual_coordinate);
			psdCluster_Z[nLayer].push_back(actualZ);

			// Search for the next maximum energy release bar
			eMaxBarPsd[nLayer] = 0; // Reset maximum energy release

			for (unsigned int barIdx = 0; barIdx < layerBarEnergyPsd[nLayer].size(); ++barIdx)
				// Check for the next unused PSD bar with a non negligible energy release
				if (layerBarUsedPsd[nLayer][barIdx] == 0 && layerBarEnergyPsd[nLayer][barIdx] > eMaxBarPsd[nLayer])
				{
					iMaxBarPsd[nLayer] = barIdx;
					eMaxBarPsd[nLayer] = layerBarEnergyPsd[nLayer][barIdx];
				}
		}
	}

	// Get the number of clusters on both X and Y
	nPsdClusters = psdCluster_idxBeg[0].size() + psdCluster_idxBeg[1].size();
	nPsdClustersX = psdCluster_idxBeg[0].size();
	nPsdClustersY = psdCluster_idxBeg[1].size();
}

std::vector<std::vector<short>> DmpPsdContainer::getPsdClusterIdxBegin()
{
	return psdCluster_idxBeg;
}

std::vector<std::vector<double>> DmpPsdContainer::getPsdClusterZ()
{
	return psdCluster_Z;
}

std::vector<std::vector<double>> DmpPsdContainer::getPsdClusterMaxE()
{
	return psdCluster_maxE;
}

std::vector<std::vector<short>> DmpPsdContainer::getPsdClusterIdxMaxE()
{
	return psdCluster_idxMaxE;
}

std::vector<std::vector<double>> DmpPsdContainer::getPsdClusterMaxECoo()
{
	return psdCluster_maxEcoordinate;
}

std::vector<double> DmpPsdContainer::getHitZ()
{
	return hitZ;
}

std::vector<short> DmpPsdContainer::getGlobalBarID()
{
	return globalBarID;
}

const unsigned int DmpPsdContainer::getPsdNclusters()
{
	return nPsdClusters;
}

const unsigned int DmpPsdContainer::getPsdNclustersX()
{
	return nPsdClustersX;
}

const unsigned int DmpPsdContainer::getPsdNclustersY()
{
	return nPsdClustersY;
}

void DmpPsdContainer::Reset()
{
	layerBarIndexPsd 					= std::vector<std::vector<short>>(DAMPE_psd_nLayers, std::vector<short>());
	layerBarNumberPsd 					= std::vector<std::vector<short>>(DAMPE_psd_nLayers, std::vector<short>());
	layerBarEnergyPsd 					= std::vector<std::vector<double>>(DAMPE_psd_nLayers, std::vector<double>());
	layerBarUsedPsd 					= std::vector<std::vector<short>>(DAMPE_psd_nLayers, std::vector<short>());
	psdCluster_idxBeg 					= std::vector<std::vector<short>>(DAMPE_psd_nLayers, std::vector<short>());
	psdCluster_length 					= std::vector<std::vector<short>>(DAMPE_psd_nLayers, std::vector<short>());
	psdCluster_idxMaxE 					= std::vector<std::vector<short>>(DAMPE_psd_nLayers, std::vector<short>());
	psdCluster_E 						= std::vector<std::vector<double>>(DAMPE_psd_nLayers, std::vector<double>());
	psdCluster_maxE	 					= std::vector<std::vector<double>>(DAMPE_psd_nLayers, std::vector<double>());
	psdCluster_maxEcoordinate 			= std::vector<std::vector<double>>(DAMPE_psd_nLayers, std::vector<double>());
	psdCluster_coordinate 				= std::vector<std::vector<double>>(DAMPE_psd_nLayers, std::vector<double>());
	psdCluster_Z 						= std::vector<std::vector<double>>(DAMPE_psd_nLayers, std::vector<double>());
	hitZ 								= std::vector<double>();
	globalBarID 						= std::vector<short>();

	nPsdClusters 						= 0;
	nPsdClustersX 						= 0;
	nPsdClustersY 						= 0;
}