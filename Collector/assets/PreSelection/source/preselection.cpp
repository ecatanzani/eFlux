#include "chain.h"
#include "bgofiducial.h"
#include "preselection.h"

#include <iostream>
#include <memory>

#include "TFile.h"

void preselection(const in_pars &input_pars) {

    auto evtch = GetFileChain(input_pars.input_path, input_pars.verbose);
    TFile* outfile = TFile::Open(input_pars.input_path.c_str(), "RECREATE");
    if (!outfile->IsOpen()) {
        std::cerr << "\n\nError writing output file [" << outfile << "]\n\n";
        exit(100);
    }

    bgofiducial_distributions(outfile, evtch, input_pars.verbose);

}