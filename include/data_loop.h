#ifndef DATA_LOOP_H
#define DATA_LOOP_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TClonesArray.h"

#include "DmpBgoContainer.h"
#include "DmpPsdContainer.h"

// DAMPESW includes
#include "DmpChain.h"
#include "DmpEvtHeader.h"
#include "DmpEvtSimuPrimaries.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpStkSiCluster.h"
#include "DmpStkTrack.h"
#include "DmpEvtPsdHits.h"

struct data_statistics
{
	unsigned int event_counter = 0;
	unsigned int events_in_range = 0;
	unsigned int events_out_range = 0;
	unsigned int events_in_saa = 0;
	unsigned int triggered_events = 0;
	unsigned int selected_events = 0;
};

extern std::vector<TH1D> evLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool verbose,
	const std::string wd);

extern bool filter_this_data_event(
	event_filter &filter,
	const std::shared_ptr<DmpEvtBgoRec> &bgorec,
	const std::shared_ptr<DmpEvtBgoHits> &bgohits,
	const cuts_conf &flux_cuts,
	const double bgoTotalE,
	DmpBgoContainer &bgoVault,
	DmpPsdContainer &psdVault,
	best_track &event_best_track,
	psd_cluster_match &clu_matching,
	psd_charge &extracted_psd_charge,
	stk_charge &extracted_stk_charge,
	const std::shared_ptr<TClonesArray> &stkclusters,
	const std::shared_ptr<TClonesArray> &stktracks,
	const data_active_cuts &active_cuts);

extern void readInputTree(
	const std::string inputPath,
	std::vector<float> &dataValues,
	TTree &dTree);

extern std::shared_ptr<TTree> getDataTree(
	const std::string inputPath,
	std::vector<float> &dataValues);

extern void branchTree(
	TTree &myDataTree,
	std::vector<float> &dataValues);

#endif