#include "chain.h"
#include "config.h"
#include "bgofiducial.h"
#include "preselection.h"

#include <iostream>
#include <memory>

#include "TFile.h"

void preselection(const in_pars &input_pars) {

    auto evtch = GetFileChain(input_pars.input_path, input_pars.verbose);
    std::shared_ptr<config> evt_config = std::make_shared<config>(input_pars.config_wd);

    bgofiducial_distributions(input_pars.output_path, input_pars.logs_dir, evtch, evt_config, input_pars.mc, input_pars.verbose);

}