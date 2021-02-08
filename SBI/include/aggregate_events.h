#ifndef AGGREGATE_EVENTS_H
#define AGGREGATE_EVENTS_H

#include <memory>
#include <string>
#include <fstream>
#include <sstream>

#include "DmpIOSvc.h"
#include "DmpCore.h"
#include "DmpFilterOrbit.h"

#include "DmpChain.h"

class event_collector
{
public:
	event_collector(
		const std::string listInputPath,
		const bool verbose = true)
	{
		input_list = listInputPath;
		verbosity = verbose;
		use_chain();
	}
	~event_collector(){};
	std::shared_ptr<DmpChain> GetChain();
	bool GetChainStatus();

private:
	std::string parse_input_file();
	void use_chain();

	std::string input_list;
	std::string tree_name = "CollectionTree";
	bool verbosity;
	std::shared_ptr<DmpChain> evt_chain;
};

#endif