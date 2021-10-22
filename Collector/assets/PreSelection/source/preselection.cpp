#include "chain.h"
#include "bgo.h"
#include "preselection.h"

#include <iostream>
#include <memory>

#include "TFile.h"

void preselection(const in_pars &input_pars) {

    auto evtch = GetFileChain(input_pars.input_path, input_pars.verbose);

    bgofiducial_distributions(input_pars, evtch);

    //lateral_and_showering_distributions(input_pars.output_wd, input_pars.logs_dir, evtch, evt_config, input_pars.mc, input_pars.verbose);

}