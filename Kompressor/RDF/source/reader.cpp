#include "reader.h"

void reader(
    const std::string wd,
    const std::string inputList,
    const std::string outputPath,
    const bool _VERBOSE,
    const bool mc)
{
    std::unique_ptr<parser> evt_parser = std::make_unique<parser>(inputList, mc, _VERBOSE);
    std::shared_ptr<config> _config = std::make_shared<config>(wd, mc);
    const double _entries = evt_parser->GetEvtTree()->GetEntries();
    if (_VERBOSE)
    {
        _config->PrintActiveFilters();
        std::cout << "Total number of events: " << _entries;
    }
    if (mc)
        mc_reader(
            evt_parser->GetEvtTree(), 
            _config, 
            _entries, 
            outputPath, 
            _VERBOSE);
}