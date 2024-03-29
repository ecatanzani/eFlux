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

    // Build acceptance histos
    auto h_gen = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Histo1D({"h_gen", "generated events; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_geometric = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                            .Filter("evtfilter_geometric_before_trigger==true")
                                            .Histo1D({"h_geometric", "generated events - geometric cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_geometric_trigger = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_geometric_before_trigger==true")
                                                .Filter("evtfilter_evt_triggered==true")
                                                .Histo1D({"h_geometric_trigger", "generated events - geometric cut + trigger; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");   

    auto h_bgo_fiducial = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_BGO_fiducial==true")
                                                .Histo1D({"h_bgo_fiducial", "generated events - BGO fiducial cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_bgo_fiducial_het = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_BGO_fiducial_HET==true")
                                                .Histo1D({"h_bgo_fiducial_het", "generated events - BGO fiducial cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_nBarLayer13 = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_nBarLayer13_cut==true")
                                                .Histo1D({"h_nBarLayer13", "generated events - nBarLayer13 cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_maxrms = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_maxRms_cut==true")
                                                .Histo1D({"h_maxrms", "generated events - max RMS cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_sumrms_low_energy = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_sumrms_low_energy_cut==true")
                                                .Histo1D({"h_sumrms_low_energy", "generated events - sumRms low energy cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_trackselection = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_track_selection_cut==true")
                                                .Histo1D({"h_trackselection", "generated events - Track Selection cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_stk_1_rm = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_stk_1rm_cut==true")
                                                .Histo1D({"h_stk_1_rm", "generated events - STK 1 RM cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_trackselection_no_3hit_recover = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_track_selection_cut_no_3hit_recover==true")
                                                .Histo1D({"h_trackselection_no_3hit_recover", "generated events - Track Selection cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_bgo_stk_selection_inside_psd_fvolume = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==true")
                                                .Histo1D({"h_bgo_stk_selection_inside_psd_fvolume", "generated events - BGO + Track Selection cut - events within PSD fiducial volume; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_bgo_stk_selection_outside_psd_fvolume = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_track_selection_cut==true")
                                                .Filter("evtfilter_psd_fiducial_volume==false")
                                                .Histo1D({"h_bgo_stk_selection_outside_psd_fvolume", "generated events - BGO + Track Selection cut - events outside PSD fiducial volume; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_psdstkmatch = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_psd_stk_match_cut==true")
                                                .Histo1D({"h_psdstkmatch", "generated events - PSD/STK match cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_psdcharge_no_one_view_recover = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_psd_charge_cut_no_single_view_recover==true")
                                                .Histo1D({"h_psdcharge_no_one_view_recover", "generated events - PSD Charge cut - no one view only recover; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");                                            
    
    auto h_psdcharge = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_psd_charge_cut==true")
                                                .Histo1D({"h_psdcharge", "generated events - PSD Charge cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_stkcharge = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_stk_charge_cut==true")
                                                .Histo1D({"h_stkcharge", "generated events - STK Charge cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    auto h_all_cut = _data_fr_selected.Define("simu_energy_gev", "simu_energy * 0.001")
                                                .Filter("evtfilter_all_cut==true")
                                                .Histo1D({"h_all_cut", "generated events - all cut; Real energy [GeV]; counts", energy_nbins, energy_binning.data()}, "simu_energy_gev");

    h_gen                                       ->Sumw2();
    h_geometric                                 ->Sumw2();
    h_geometric_trigger                         ->Sumw2();
    h_bgo_fiducial                              ->Sumw2();
    h_bgo_fiducial_het                          ->Sumw2();
    h_nBarLayer13                               ->Sumw2();
    h_maxrms                                    ->Sumw2();
    h_sumrms_low_energy                         ->Sumw2();
    h_trackselection                            ->Sumw2();
    h_trackselection_no_3hit_recover            ->Sumw2();
    h_stk_1_rm                                  ->Sumw2();
    h_bgo_stk_selection_inside_psd_fvolume      ->Sumw2();
    h_bgo_stk_selection_outside_psd_fvolume     ->Sumw2();
    h_psdstkmatch                               ->Sumw2();
    h_psdcharge_no_one_view_recover             ->Sumw2();
    h_psdcharge                                 ->Sumw2();
    h_stkcharge                                 ->Sumw2();  
    h_all_cut                                   ->Sumw2();

    TFile* outfile = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "Error writing output acceptance file: [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    h_gen                                       ->Write();
    h_geometric                                 ->Write();
    h_geometric_trigger                         ->Write();
    h_bgo_fiducial                              ->Write();
    h_bgo_fiducial_het                          ->Write();
    h_nBarLayer13                               ->Write();
    h_maxrms                                    ->Write();
    h_sumrms_low_energy                         ->Write();  
    h_trackselection                            ->Write();
    h_trackselection_no_3hit_recover            ->Write();
    h_stk_1_rm                                  ->Write();
    h_bgo_stk_selection_inside_psd_fvolume      ->Write();
    h_bgo_stk_selection_outside_psd_fvolume     ->Write();
    h_psdstkmatch                               ->Write();
    h_psdcharge_no_one_view_recover             ->Write();
    h_psdcharge                                 ->Write();
    h_stkcharge                                 ->Write();
    h_all_cut                                   ->Write();

    outfile->Close();

    if (input_args.verbose)
        std::cout << "Output file has been written: [" << input_args.output_path << "]\n\n";
}