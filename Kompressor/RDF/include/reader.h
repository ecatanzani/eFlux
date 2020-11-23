#ifndef READER_H
#define READER_H

#include <memory>
#include <string>

#include "TChain.h" 

#include "config.h"
#include "list_parser.h"

#define nthreads 20

void reader(
    const std::string wd,
    const std::string inputList,
    const std::string outputPath,
    const bool _VERBOSE,
    const bool mc);

void mc_reader(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE);

#endif
