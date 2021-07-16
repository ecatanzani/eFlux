#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <memory>

#include "main.h"
#include "parser.h"
#include "energy_config.h"

extern void buildEfficiency(const in_args input_args);
extern void buildMCErriciency(
    const in_args input_args,
    std::shared_ptr<parser> evt_parser, 
    std::shared_ptr<energy_config> mc_energy_config);
extern void buildDATAErriciency(
    const in_args input_args,
    std::shared_ptr<parser> evt_parser, 
    std::shared_ptr<energy_config> mc_energy_config);

#endif