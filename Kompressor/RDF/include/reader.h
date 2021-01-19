#ifndef READER_H
#define READER_H

#include <memory>
#include <string>
#include <memory>

#include "TChain.h" 

#include "main.h"
#include "config.h"
#include "list_parser.h"

extern void reader(in_args input_args);

extern void mc_reader(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    const std::string fit_tree_path,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE);

#endif
