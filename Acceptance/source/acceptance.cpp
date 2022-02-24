#include "config.h"
#include "parser.h"
#include "acceptance.h"
#include "energy_config.h"

#include <memory>

#include "TFile.h"
#include <ROOT/RDataFrame.hxx>

#define _NO_STK_CHARGE true

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
    
    // Filter events in the interesting energy range
    auto energyFilter = [&mc_energy_config](const double energy) -> bool 
    {
        auto status = false;
        auto _gev = 0.001;
        if (energy*_gev >= mc_energy_config->GetMinEvtEnergy() && energy*_gev <= mc_energy_config->GetMaxEvtEnergy())
            status = true;
        return status; 
    };

    auto _data_fr_selected = _data_fr.Filter(energyFilter, {"simu_energy"}).Filter("evtfilter_correct_bgo_reco==1");

    // Build acceptance histos
    auto h_gen = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Histo1D({"h_gen", "generated events; Real energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_geometric = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                            .Filter("evtfilter_geometric_before_trigger==true")
                                            .Histo1D({"h_geometric", "generated events - geometric cut; Real energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_geometric_trigger = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_geometric_before_trigger==true")
                                                .Filter("evtfilter_evt_triggered==true")
                                                .Histo1D({"h_geometric_trigger", "generated events - geometric cut + trigger; Real energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");                                        
    auto h_bgo_fiducial = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_BGO_fiducial==true")
                                                .Histo1D({"h_bgo_fiducial", "generated events - BGO fiducial cut; Real energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    
    auto h_all_cut = _NO_STK_CHARGE ? _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_psd_charge_cut==true")
                                                .Histo1D({"h_all_cut", "generated events - all cut; Real energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev") :  _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_all_cut==true")
                                                .Histo1D({"h_all_cut", "generated events - all cut; Real energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");

    h_gen->Sumw2();
    h_geometric->Sumw2();
    h_geometric_trigger->Sumw2();
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
    h_geometric_trigger->Write();
    h_bgo_fiducial->Write();
    h_all_cut->Write();

    outfile->Close();

    if (input_args.verbose)
        std::cout << "Output file has been written: [" << input_args.output_path << "]\n\n";
}