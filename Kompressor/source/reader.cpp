#include "main.h"
#include "utils.h"
#include "config.h"
#include "mc_tuple.h"
#include "mc_histos.h"
#include "environment.h"
#include "list_parser.h"

void reader(
    const std::string wd,
    const std::string inputList,
    const std::string outputPath,
    const bool _VERBOSE,
    const bool mc)
{
    std::unique_ptr<parser> evt_parser = std::make_unique<parser>(inputList, mc, _VERBOSE);
    std::shared_ptr<config> _config = std::make_shared<config>(wd, mc);
    if (_VERBOSE)
	{
		_config->PrintActiveFilters();
		std::cout << "Reading events..." << std::endl;
	}
    if (mc)
        mc_reader(evt_parser->GetEvtTree(), _config, outputPath, _VERBOSE);
}

void mc_reader(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    const std::string outputPath,
    const bool _VERBOSE)
{
#if _PARALLER
    ROOT::EnableImplicitMT(nthreads);
#endif
    int kStep = 10;
    std::unique_ptr<mc_tuple> _tuple = std::make_unique<mc_tuple>(evtch);
    _tuple->InitHistos(_config->GetEnergyBinning());
    auto _entries = evtch->GetEntries();
    for (auto i : ROOT::TSeqUL(_entries)) 
    {
        if (_VERBOSE)
            UpdateProcessStatus(i, kStep, _entries);
        _tuple->GetEntry(i);
    }    
    _tuple->WriteHistos(outputPath);
}