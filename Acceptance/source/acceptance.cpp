#include "config.h"
#include "parser.h"
#include "acceptance.h"
#include "energy_config.h"

#include <memory>

#include "TFile.h"
#include <ROOT/RDataFrame.hxx>

void buildAcceptance(const in_args input_args)
{
    std::shared_ptr<parser> evt_parser = std::make_shared<parser>(input_args.input_list, input_args.verbose);
    std::shared_ptr<config> mc_config = std::make_shared<config>(input_args.wd);
    std::shared_ptr<energy_config> mc_energy_config = std::make_shared<energy_config>(input_args.wd);

    if (input_args.verbose)
    {
        mc_config->PrintActiveFilters();
        std::cout << "\n\nAnalysis running...\n\n";
    }

    auto energy_binning = mc_energy_config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;

    // Build RDF
    ROOT::EnableImplicitMT(input_args.threads);
    ROOT::RDataFrame _data_fr(*(evt_parser->GetEvtTree()));

    // Build acceptance histos
    auto h_gen = _data_fr.Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Histo1D({"h_gen", "generated events; Real energy [GeV]", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_geometric = _data_fr.Define("simu_energy_gev", "simu_energy * 0.001")
                                            .Filter("evtfilter_geometric==true")
                                            .Histo1D({"h_geometric", "generated events - geometric cut; Real energy [GeV]", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_bgo_fiducial = _data_fr.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_BGO_fiducial==true")
                                                .Histo1D({"h_bgo_fiducial", "generated events - BGO fiducial cut; Real energy [GeV]", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_all_cut = _data_fr.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_all_cut==true")
                                                .Histo1D({"h_all_cut", "generated events - all cut; Real energy [GeV]", energy_nbins, &energy_binning[0]}, "simu_energy_gev");

    h_gen->Sumw2();
    h_geometric->Sumw2();
    h_bgo_fiducial->Sumw2();
    h_all_cut->Sumw2();

    TFile* outfile = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "Error writing output acceptance file: [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    h_gen->Write();
    h_geometric->Write();
    h_bgo_fiducial->Write();
    h_all_cut->Write();

    outfile->Close();

    if (input_args.verbose)
        std::cout << "Output file has been written: [" << input_args.output_path << "]\n\n";
}