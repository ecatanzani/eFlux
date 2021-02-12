#ifndef AGGREGATE_EVENTS_H
#define AGGREGATE_EVENTS_H

#include <memory>
#include <string>
#include <fstream>
#include <sstream>

#include "DmpIOSvc.h"
#include "DmpCore.h"
#include "DmpFilterOrbit.h"

#include "TChain.h"

class event_collector
{
public:
	event_collector(
		const std::string listInputPath,
		const bool verbose = true,
		const bool mc = false)
	{
		input_list = listInputPath;
		verbosity = verbose;
		simu_evt = mc;
		use_chain();
	}
	~event_collector(){};
	std::shared_ptr<TChain> GetChain();
	bool GetChainStatus();

private:
	std::string parse_input_file();
	void use_chain();

	std::string input_list;
	std::string tree_name = "CollectionTree";
	bool verbosity;
	bool simu_evt;
	std::shared_ptr<TChain> evt_chain;
};

extern std::shared_ptr<TChain> aggregateEventsTChain(
	const std::string listInputPath,
	const bool verbose);

extern std::shared_ptr<TChain> aggregateDataEventsTChain(
	const std::string listInputPath,
	const bool verbose,
	const bool skimmed);

extern std::shared_ptr<TChain> aggregateTupleDataEventsTChain(
	const std::string listInputPath,
	std::string &year,
	std::string &month,
	const bool verbose);

#endif