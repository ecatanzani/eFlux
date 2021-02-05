#include "aggregate_events.h"

void event_collector::use_chain()
{
	evt_chain = std::make_shared<DmpChain>(tree_name.c_str());
	std::istringstream input_stream(parse_input_file());
	std::string tmp_str;
	while (input_stream >> tmp_str)
	{
		if (!simu_evt)
		{
			gIOSvc->Set("InData/Read", tmp_str.c_str());
			if (_f_file)
			{
				auto yIdx = tmp_str.find("/2A/") + 4;
				auto mIdx = yIdx + 4;
				data_year = tmp_str.substr(yIdx, 4);
				data_month = tmp_str.substr(mIdx, 2);
				_f_file = false;
			}
		}
		evt_chain->Add(tmp_str.c_str());
		if (verbosity)
			std::cout << "\nAdding " << tmp_str << " to the chain ...";
	}
}

std::string event_collector::parse_input_file()
{
	std::ifstream input_file(input_list.c_str());
	if (!input_file.is_open())
	{
		std::cerr << "\n\nError (100) reading input file list...[" << input_list << "]" << std::endl;
		exit(100);
	}
	std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

std::shared_ptr<TChain> event_collector::GetChain()
{
	if (evt_chain)
		return evt_chain;
	else
		return std::shared_ptr<TChain>(nullptr);
}

bool event_collector::GetChainStatus()
{
	if (evt_chain)
		return true;
	else
		return false;
}