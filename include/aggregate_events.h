#ifndef AGGREGATE_EVENTS_H
#define AGGREGATE_EVENTS_H

#include <memory>
#include <string>
#include <fstream>
#include <sstream>

#include "DmpChain.h"

extern std::shared_ptr<DmpChain> aggregateEventsDmpChain(
    const std::string listInputPath,
    const bool verbose);

extern std::shared_ptr<TChain> aggregateEventsTChain(
    const std::string listInputPath,
    const bool verbose);

#endif