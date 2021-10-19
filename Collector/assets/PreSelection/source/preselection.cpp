#include "chain.h"
#include "bgofiducial.h"
#include "preselection.h"

#include <iostream>
#include <memory>

#include "TFile.h"

void preselection(const in_pars &input_pars) {

    auto evtch = GetFileChain(input_pars.input_path, input_pars.verbose);

    bgofiducial_distributions(input_pars.output_path, evtch, input_pars.verbose);

}