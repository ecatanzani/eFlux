#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include <string>
#include <vector>
#include <memory>

#include "TF1.h"
#include "TFile.h"
#include "TClonesArray.h"

#include "anyoption.h"
#include "data_cuts.h"

#include "DmpEvtSimuPrimaries.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtBgoHits.h"
#include "DmpStkSiCluster.h"
#include "DmpStkTrack.h"

#include "DmpBgoContainer.h"
#include "DmpPsdContainer.h"

struct mc_ancillary_cuts
{
    bool compute_proton_background = false;
};

struct mc_statistics
{
    unsigned int event_counter = 0;
    unsigned int generated_events_in_range = 0;
    unsigned int generated_events_out_range = 0;
    unsigned int triggered_events = 0;
    unsigned int selected_events = 0;
};

extern void computeAcceptance(
    const std::string accInputPath,
    const bool verbose,
    const bool pedantic,
    const std::string outputPath,
    AnyOption &opt,
    const std::string wd);

extern void buildAcceptance(
    const std::string accInputPath,
    const bool verbose,
    TFile &outFile,
    const std::string wd);

extern void buildAcceptance_vector(
    const std::string accInputPath,
    const bool verbose,
    TFile &outFile,
    const std::string wd);

extern TF1 readAcceptance(
    TFile &outFile,
    const bool verbose,
    const std::string accInputPath);

extern void computeAcceptance(
    const std::string accInputPath,
    const bool verbose,
    const bool pedantic);

extern std::vector<float> load_acceptance_struct(
    cuts_conf &acceptance_cuts,
    data_active_cuts &active_cuts,
    mc_ancillary_cuts &ancillary_cuts,
    const std::string wd);

extern void print_filter_status(data_active_cuts active_cuts);

extern bool filter_this_mc_event(
    event_filter &filter,
    const std::shared_ptr<DmpEvtSimuPrimaries> &simu_primaries,
    const std::shared_ptr<DmpEvtBgoRec> &bgorec,
    const std::shared_ptr<DmpEvtBgoHits> &bgohits,
    const cuts_conf &acceptance_cuts,
    const double bgoTotalE,
	DmpBgoContainer &bgoVault,
	DmpPsdContainer &psdVault,
    psd_charge &extracted_psd_charge,
	stk_charge &extracted_stk_charge,
    const std::shared_ptr<TClonesArray> &stkclusters,
    const std::shared_ptr<TClonesArray> &stktracks,
	const data_active_cuts &active_cuts,
    const mc_ancillary_cuts &ancillary_cuts);

#endif