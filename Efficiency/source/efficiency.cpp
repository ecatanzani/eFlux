#include "config.h"
#include "parser.h"
#include "efficiency.h"
#include "energy_config.h"

#include <memory>

#include "TFile.h"
#include <ROOT/RDataFrame.hxx>

void buildEfficiency(const in_args input_args)
{
    std::shared_ptr<parser> evt_parser = std::make_shared<parser>(input_args.input_list, input_args.verbose);
    std::shared_ptr<config> mc_config = std::make_shared<config>(input_args.wd);
    std::shared_ptr<energy_config> mc_energy_config = std::make_shared<energy_config>(input_args.wd);

    if (input_args.verbose)
    {
        mc_config->PrintActiveFilters();
        std::cout << "\n\nAnalysis running...\n\n";
    }

    if (input_args.mc) buildMCErriciency(input_args, evt_parser, mc_energy_config);
    else buildDATAErriciency(input_args, evt_parser, mc_energy_config);
}

void buildMCErriciency(
    const in_args input_args,
    std::shared_ptr<parser> evt_parser, 
    std::shared_ptr<energy_config> mc_energy_config)
{
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

    auto _data_fr_selected = _data_fr.Filter(energyFilter, {"simu_energy"});
    
    // Build histos
    auto h_geometric = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                        .Filter("evtfilter_geometric_before_trigger==true")
                                        .Histo1D({"h_geometric", "Geometric factor; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_geometric_trigger = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_geometric_before_trigger==true")
                                                .Filter("evtfilter_evt_triggered==true")
                                                .Histo1D({"h_geometric_trigger", "Geometric factor + trigger; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_trigger = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                        .Filter("evtfilter_evt_triggered==true")
                                        .Histo1D({"h_trigger", "Trigger", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_maxElayer = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                        .Filter("evtfilter_BGO_fiducial_maxElayer_cut==true")
                                        .Histo1D({"h_maxElayer", "maxElayer cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_maxBarlayer = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                            .Filter("evtfilter_BGO_fiducial_maxBarLayer_cut==true")
                                            .Histo1D({"h_maxBarlayer", "maxBarLayer cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_BGOTrackContainment = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                    .Filter("evtfilter_BGO_fiducial_BGOTrackContainment_cut==true")
                                                    .Histo1D({"h_BGOTrackContainment", "BGO track containment cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_bgo_fiducial = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                            .Filter("evtfilter_BGO_fiducial==true")
                                            .Histo1D({"h_bgo_fiducial", "BGO fiducial; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_nbarlayer13 = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                            .Filter("evtfilter_nBarLayer13_cut==true")
                                            .Histo1D({"h_nbarlayer13", "nBar layer 13 cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_maxrms = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                        .Filter("evtfilter_maxRms_cut==true")
                                        .Histo1D({"h_maxrms", "max RMS cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_track_selection = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_track_selection_cut==true")
                                                .Histo1D({"h_track_selection", "Track selection cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_psd_stk_match = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_psd_stk_match_cut==true")
                                                .Histo1D({"h_psd_stk_match", "PSD-STK match cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_psd_charge = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                            .Filter("evtfilter_psd_charge_cut==true")
                                            .Histo1D({"h_psd_charge", "PSD charge cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_stk_charge = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                            .Filter("evtfilter_stk_charge_cut==true")
                                            .Histo1D({"h_stk_charge", "STK charge cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");
    auto h_all_cuts = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                        .Filter("evtfilter_all_cut==true")
                                        .Histo1D({"h_all_cuts", "All cuts; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "simu_energy_gev");

    h_geometric->Sumw2();
    h_geometric_trigger->Sumw2();
    h_trigger->Sumw2();
    h_maxElayer->Sumw2();
    h_maxBarlayer->Sumw2();
    h_BGOTrackContainment->Sumw2();
    h_bgo_fiducial->Sumw2();
    h_nbarlayer13->Sumw2();
    h_maxrms->Sumw2();
    h_track_selection->Sumw2();
    h_psd_stk_match->Sumw2();
    h_psd_charge->Sumw2();
    h_stk_charge->Sumw2();
    h_all_cuts->Sumw2();

    TFile* outfile = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "Error writing output acceptance file: [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    h_geometric->Write();
    h_geometric_trigger->Write();
    h_trigger->Write();
    h_maxElayer->Write();
    h_maxBarlayer->Write();
    h_BGOTrackContainment->Write();
    h_bgo_fiducial->Write();
    h_nbarlayer13->Write();
    h_maxrms->Write();
    h_track_selection->Write();
    h_psd_stk_match->Write();
    h_psd_charge->Write();
    h_stk_charge->Write();
    h_all_cuts->Write();

    outfile->Close();

    if (input_args.verbose)
        std::cout << "Output file has been written: [" << input_args.output_path << "]\n\n";
}

void buildDATAErriciency(
    const in_args input_args,
    std::shared_ptr<parser> evt_parser, 
    std::shared_ptr<energy_config> mc_energy_config)
{
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

    auto _data_fr_selected = _data_fr.Filter(energyFilter, {"energy_corr"});

    // Build histos
    auto h_trigger = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                        .Filter("evtfilter_evt_triggered==true")
                                        .Histo1D({"h_trigger", "Trigger", energy_nbins, &energy_binning[0]}, "energy_corr_gev");
    auto h_maxElayer = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                            .Filter("evtfilter_BGO_fiducial_maxElayer_cut==true")
                                            .Histo1D({"h_maxElayer", "maxElayer cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "energy_corr_gev");
    auto h_maxBarlayer = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                            .Filter("evtfilter_BGO_fiducial_maxBarLayer_cut==true")
                                            .Histo1D({"h_maxBarlayer", "maxBarLayer cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "energy_corr_gev");
    auto h_BGOTrackContainment = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                                    .Filter("evtfilter_BGO_fiducial_BGOTrackContainment_cut==true")
                                                    .Histo1D({"h_BGOTrackContainment", "BGO track containment cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "energy_corr_gev");
    auto h_bgo_fiducial = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                            .Filter("evtfilter_BGO_fiducial==true")
                                            .Histo1D({"h_bgo_fiducial", "BGO fiducial; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "energy_corr_gev");
    auto h_nbarlayer13 = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                            .Filter("evtfilter_nBarLayer13_cut==true")
                                            .Histo1D({"h_nbarlayer13", "nBar layer 13 cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "energy_corr_gev");
    auto h_maxrms = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                        .Filter("evtfilter_maxRms_cut==true")
                                        .Histo1D({"h_maxrms", "max RMS cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "energy_corr_gev");
    auto h_track_selection = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                                .Filter("evtfilter_track_selection_cut==true")
                                                .Histo1D({"h_track_selection", "Track selection cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "energy_corr_gev");
    auto h_psd_stk_match = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                            .Filter("evtfilter_psd_stk_match_cut==true")
                                            .Histo1D({"h_psd_stk_match", "PSD-STK match cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "energy_corr_gev");
    auto h_psd_charge = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                            .Filter("evtfilter_psd_charge_cut==true")
                                            .Histo1D({"h_psd_charge", "PSD charge cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "energy_corr_gev");
    auto h_stk_charge = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                            .Filter("evtfilter_stk_charge_cut==true")
                                            .Histo1D({"h_stk_charge", "STK charge cut; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "energy_corr_gev");
    auto h_all_cuts = _data_fr_selected.Define("energy_corr_gev", "energy_corr * 0.001")
                                        .Filter("evtfilter_all_cut==true")
                                        .Histo1D({"h_all_cuts", "All cuts; Real Energy [GeV]; counts", energy_nbins, &energy_binning[0]}, "energy_corr_gev");

    h_trigger->Sumw2();
    h_maxElayer->Sumw2();
    h_maxBarlayer->Sumw2();
    h_BGOTrackContainment->Sumw2();
    h_bgo_fiducial->Sumw2();
    h_nbarlayer13->Sumw2();
    h_maxrms->Sumw2();
    h_track_selection->Sumw2();
    h_psd_stk_match->Sumw2();
    h_psd_charge->Sumw2();
    h_stk_charge->Sumw2();
    h_all_cuts->Sumw2();

    TFile* outfile = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "Error writing output acceptance file: [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    h_trigger->Write();
    h_maxElayer->Write();
    h_maxBarlayer->Write();
    h_BGOTrackContainment->Write();
    h_bgo_fiducial->Write();
    h_nbarlayer13->Write();
    h_maxrms->Write();
    h_track_selection->Write();
    h_psd_stk_match->Write();
    h_psd_charge->Write();
    h_stk_charge->Write();
    h_all_cuts->Write();

    outfile->Close();

    if (input_args.verbose)
        std::cout << "Output file has been written: [" << input_args.output_path << "]\n\n";
}