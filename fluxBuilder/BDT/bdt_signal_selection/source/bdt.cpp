#include "bdt.h"
#include "list_parser.h"
#include "energy_config.h"

#include <memory>

#include "TFile.h"
#include <ROOT/RDataFrame.hxx>

void bdt_selection(in_args input_args) {

    std::shared_ptr<energy_config> config = std::make_shared<energy_config>(input_args.energy_config_file);
    std::shared_ptr<parser> evt_parser = std::make_unique<parser>(input_args.input_list, input_args.verbose);

    auto energy_binning = config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;

    ROOT::EnableImplicitMT(input_args.threads);
    ROOT::RDataFrame _data_fr(*evt_parser->GetEvtTree());

    std::cout << "\n\n**** Filter statistics ****\n";
    std::cout << "***************************\n";
    std::cout << "\nTotal events: " << *(_data_fr.Count());
    std::cout << "\n\n***************************";

    if (input_args.verbose) std::cout << "\n\nAnlysis running...\n\n";

    auto h_bdt_selection = _data_fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Histo1D<double>({"h_bdt_selection", "BDT selected events; BGO Corr energy [GeV]; entries", energy_nbins, &energy_binning[0]}, "corr_energy_gev");

    TFile* outfile = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (!outfile->IsOpen()) {
        std::cerr << "\n\nError writing output file... [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    h_bdt_selection->Write();

    outfile->Close();
}