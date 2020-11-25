#include "reader.h"
#include "binning.h"
#include "DAMPE_geo_structure.h"

#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TVector3.h"
#include "TDirectory.h"

void reader(
    const std::string wd,
    const std::string inputList,
    const std::string outputPath,
    const bool _VERBOSE,
    const bool mc)
{
    std::unique_ptr<parser> evt_parser = std::make_unique<parser>(inputList, mc, _VERBOSE);
    std::shared_ptr<config> _config = std::make_shared<config>(wd, mc);
    const double _entries = evt_parser->GetEvtTree()->GetEntries();
    if (_VERBOSE)
    {
        _config->PrintActiveFilters();
        std::cout << "Total number of events: " << _entries;
    }
    if (mc)
        mc_reader(evt_parser->GetEvtTree(), _config, _entries, outputPath, _VERBOSE);
}

void mc_reader(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE)
{
    // Enable multithreading
    ROOT::EnableImplicitMT();
    // Create RDF
    ROOT::RDataFrame _data_fr(*evtch);
    // Extract the energy binning
    auto energy_binning = _config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;
    double _gev = 0.001;
    // Create the RDF with the corrected energy weight
    auto GetEnergyBin = [=](double energy) -> int { 
        int bin_idx=0;
        for (; bin_idx<energy_nbins; ++bin_idx)
            if (energy * _gev >= energy_binning[bin_idx] && energy * _gev < energy_binning[bin_idx+1])
                break;
        return bin_idx+1; };
    auto _fr_weight_patch = _data_fr.Define("simu_energy_w_shift", [&energy_binning] { return pow(energy_binning[0], 2); })
                                .Define("simu_energy_w_pathc", "simu_energy_w_shift*simu_energy_w")
                                .Define("energy_bin", GetEnergyBin, {"energy_corr"});
    // Create the RDF for BGO analysis
    auto _fr_bgo_analysis = _fr_weight_patch.Filter("evtfilter_good_event==true");
    // Create the RDF for the preselected events
    auto _fr_preselected = _fr_bgo_analysis.Filter("evtfilter_all_cut==true");
    // Create the RDF for STK analysis
    auto _fr_stk_analysis = _fr_bgo_analysis.Filter("evtfilter_track_selection_cut==true");
    // Create the RDF for PSD charge analysis
    auto _fr_psd_charge_analysis = _fr_stk_analysis.Filter("evtfilter_psd_charge_measurement==true");
    // Create the RDF for STK charge analysis
    auto _fr_stk_charge_analysis = _fr_stk_analysis.Filter("evtfilter_stk_charge_measurement==true");

    std::cout << "\nBGO analysis events: " << *(_fr_bgo_analysis.Count());
    std::cout << "\nSTK analysis events: " << *(_fr_stk_analysis.Count());
    std::cout << "\nPSD charge analysis events: " << *(_fr_psd_charge_analysis.Count());
    std::cout << "\nSTK charge analysis events: " << *(_fr_stk_charge_analysis.Count());
    std::cout << "\nPreselected events: " << *(_fr_preselected.Count());

    if (_VERBOSE)
        std::cout << "\n\nAnlysis running..." << std::endl;

    // Extract BGO histos
    auto sumRms_bins = createLogBinning(10, 2e+3, 1e+2);
    auto xtrl_bins = createLinearBinning(0, 150, 1e+2);
    auto cosine_bins = createLinearBinning(0, 1, 1e+2);
    auto energy_ratio_bins = createLinearBinning(-1, 1, 1e+2);
    auto flast_binning = createLogBinning(1e-5, 2e-1, 1e+2);

    auto h_BGOrec_raw_energy = _fr_bgo_analysis.Define("raw_energy_gev", "energy * 0.001")
                                   .Histo1D<double, double>({"h_BGOrec_raw_energy", "BGO Raw energy", energy_nbins, &energy_binning[0]}, "raw_energy_gev", "simu_energy_w_pathc");
    auto h_BGOrec_corr_energy = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Histo1D<double, double>({"h_BGOrec_corr_energy", "BGO Corr energy", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "simu_energy_w_pathc");
    auto h_BGOrec_cosine = _fr_bgo_analysis.Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                               .Histo1D<double, double>({"h_BGOrec_cosine", "BGO Reco cosine", 100, 0, 1}, "bgorec_cosine", "simu_energy_w_pathc");
    auto h_BGOrec_cosine_preseltion = _fr_preselected.Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                          .Histo1D<double, double>({"h_BGOrec_cosine_preseltion", "BGO Reco cosine - Preselection", 100, 0, 1}, "bgorec_cosine", "simu_energy_w_pathc");

    auto h_BGOrec_sumRms_bin_cosine2D_20_100 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                   .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                                   .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                   .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_20_100", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_bin_cosine2D_100_250 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                    .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                                    .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                    .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_100_250", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_bin_cosine2D_250_500 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                    .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                                    .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                    .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_250_500", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_bin_cosine2D_500_1000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                     .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                                     .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                     .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_500_1000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_bin_cosine2D_1000_3000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                      .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                                      .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                      .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_1000_3000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_bin_cosine2D_3000_10000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                       .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                                       .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                       .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_3000_10000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_bin_cosine2D_10000_20000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                                        .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                        .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_10000_20000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");

    auto h_BGOrec_sumRms_flast_20_100 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                            .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_20_100", "sumRms vs F_{last} correlation - 20 GeV - 100 GeV; sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_100_250 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                             .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                             .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_100_250", "sumRms vs F_{last} correlation - 100 GeV - 250 GeV;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_250_500 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                             .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                             .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_250_500", "sumRms vs F_{last} correlation - 250 GeV - 500 GeV;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_500_1000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                              .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                              .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_500_1000", "sumRms vs F_{last} correlation - 500 GeV - 1 TeV;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_1000_3000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                               .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                               .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_1000_3000", "sumRms vs F_{last} correlation - 1 TeV - 3 TeV;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_3000_10000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                                .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_3000_10000", "sumRms vs F_{last} correlation - 3 TeV - 10 TeV ;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_10000_20000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                 .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                                 .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_10000_20000", "sumRms vs F_{last} correlation - 10 TeV - 20 TeV ;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");

    auto h_BGOrec_sumRms_flast_13_20_100 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                               .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                               .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_20_100", "sumRms vs F_{13} correlation - 20 GeV - 100 GeV ;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_13_100_250 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                                .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_100_250", "sumRms vs F_{13} correlation - 100 GeV - 250 GeV ;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_13_250_500 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                                .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_250_500", "sumRms vs F_{13} correlation - 250 GeV - 500 GeV ;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_13_500_1000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                 .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                                 .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_500_1000", "sumRms vs F_{13} correlation - 500 GeV - 1 TeV ;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_13_1000_3000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                  .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                                  .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_1000_3000", "sumRms vs F_{13} correlation - 1 TeV - 3 TeV ;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_13_3000_10000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                   .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                                   .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_3000_10000", "sumRms vs F_{13} correlation - 3 TeV - 10 TeV ;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_BGOrec_sumRms_flast_13_10000_20000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                    .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                                    .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_10000_20000", "sumRms vs F_{13} correlation - 10 TeV - 20 TeV ;sumRms [mm]; F_{last}", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");

    auto h_BGOrec_slopeX = _fr_bgo_analysis.Histo1D<double, double>({"h_BGOrec_slopeX", "BGOrec Slope X", 200, -10, 10}, "BGOrec_slopeX", "simu_energy_w_pathc");
    auto h_BGOrec_slopeY = _fr_bgo_analysis.Histo1D<double, double>({"h_BGOrec_slopeY", "BGOrec Slope Y", 200, -10, 10}, "BGOrec_slopeY", "simu_energy_w_pathc");
    auto h_BGOrec_interceptX = _fr_bgo_analysis.Histo1D<double, double>({"h_BGOrec_interceptX", "BGOrec Intercept X", 100, -500, 500}, "BGOrec_interceptX", "simu_energy_w_pathc");
    auto h_BGOrec_interceptY = _fr_bgo_analysis.Histo1D<double, double>({"h_BGOrec_interceptY", "BGOrec Intercept Y", 100, -500, 500}, "BGOrec_interceptY", "simu_energy_w_pathc");

    auto h_xtrl_energy_int = _fr_bgo_analysis.Histo1D<double>({"h_xtrl_energy_int", "XTRL", 200, 0, 150}, "xtrl", "simu_energy_w_pathc");
    auto h_xtrl = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                      .Histo2D<double, double, double>({"h_xtrl", "XTRL", energy_nbins, &energy_binning[0], (int)xtrl_bins.size() - 1, &xtrl_bins[0]}, "corr_energy_gev", "xtrl", "simu_energy_w_pathc");

    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_cosine_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_sumRms_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_sumRms_mean_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_sumRms_bin_weight(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_sumRms_sumRms_weight(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_sumRms_bin_cosine(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_sumRms_bin_cosine2D(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_sumRms_flast_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_sumRms_flast_13_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_max_energy_ratio(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ratio_last(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ratio_13(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ratio_last_cosine(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_last_layer(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_hits(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_xtrl_bin(energy_nbins);
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_rms_layer(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_energy_ratio_layer(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_energy_ratio_1R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_energy_ratio_2R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_energy_ratio_3R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_energy_ratio_5R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));

    std::vector<std::shared_ptr<TH2D>> h_BGOrec_shower_profile(energy_nbins);
    std::vector<std::shared_ptr<TH2D>> h_BGOrec_shower_profile_upto_09(energy_nbins);
    std::vector<std::shared_ptr<TH2D>> h_BGOrec_shower_profile_from_09(energy_nbins);
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        h_BGOrec_shower_profile[bin_idx - 1] = std::make_shared<TH2D>((std::string("h_BGOrec_shower_profile_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Shower Profile - bin ") + std::to_string(bin_idx)).c_str(), 13, 0, 13, 100, 0, 1);
        h_BGOrec_shower_profile_upto_09[bin_idx - 1] = std::make_shared<TH2D>((std::string("h_BGOrec_shower_profile_upto_09_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Shower Profile cos(#theta) < 0.9 - bin ") + std::to_string(bin_idx)).c_str(), 13, 0, 13, 100, 0, 1);
        h_BGOrec_shower_profile_from_09[bin_idx - 1] = std::make_shared<TH2D>((std::string("h_BGOrec_shower_profile_from_09_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Shower Profile cos(#theta) > 0.9 - bin ") + std::to_string(bin_idx)).c_str(), 13, 0, 13, 100, 0, 1);
    }

    auto computeSumRmsWeight = [=](std::vector<double> &energy_layer, std::vector<double> &rms_layer, double raw_energy) {
        double sumRms_weight = 0;
        for (int ly = 0; ly<DAMPE_bgo_nLayers; ++ly)
            sumRms_weight += energy_layer[ly]*rms_layer[ly];
        return sumRms_weight/raw_energy; };
    auto GetMaxEnergyRatio = [](std::vector<double> &energy_fraction) -> double { 
        auto it_max = std::max_element(energy_fraction.begin(), energy_fraction.end());
        return energy_fraction[std::distance(energy_fraction.begin(), it_max)]; };

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        auto bin_filter = [=](int energy_bin) -> bool { return energy_bin == bin_idx; };
        h_BGOrec_cosine_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                               .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                               .Histo1D<double, double>({(std::string("h_BGOrec_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Reco cosine - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "bgorec_cosine", "simu_energy_w_pathc");
        h_BGOrec_sumRms_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                               .Histo1D<double, double>({(std::string("h_BGOrec_sumRms_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 3000}, "sumRms", "simu_energy_w_pathc");
        h_BGOrec_sumRms_mean_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                    .Define("num_layers", [](int last_layer) -> int { return last_layer + 1; }, {"lastBGOLayer"})
                                                    .Define("sumRms_eff", "sumRms/(lastBGOLayer+1)")
                                                    .Histo1D<double, double>({(std::string("h_BGOrec_sumRms_mean_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms - bin ") + std::to_string(bin_idx)).c_str(), 50, 0, 100}, "sumRms_eff", "simu_energy_w_pathc");
        h_BGOrec_sumRms_bin_weight[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                      .Define("sumRms_w", computeSumRmsWeight, {"eLayer", "rmsLayer", "energy"})
                                                      .Histo1D<double, double>({(std::string("h_BGOrec_sumRms_weight_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms weighted - bin ") + std::to_string(bin_idx)).c_str(), 50, 0, 100}, "sumRms_w", "simu_energy_w_pathc");
        h_BGOrec_sumRms_sumRms_weight[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                         .Define("sumRms_w", computeSumRmsWeight, {"eLayer", "rmsLayer", "energy"})
                                                         .Define("sumRms_eff", "sumRms/(lastBGOLayer+1)")
                                                         .Histo2D<double, double, double>({(std::string("h_BGOrec_sumRms_sumRms_weight_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms_{mean} vs sumRms_{weighted} - bin ") + std::to_string(bin_idx) + std::string("; sumRms_{weighted} [mm]; sumRms_{mean} [mm]")).c_str(), 50, 0, 100, 50, 0, 100}, "sumRms_w", "sumRms_eff", "simu_energy_w_pathc");
        h_BGOrec_sumRms_bin_cosine[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                      .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                      .Define("sumRms_cosine", "sumRms/bgorec_cosine")
                                                      .Histo1D<double, double>({(std::string("h_BGOrec_sumRms_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; sumRms/cos(#theta); counts")).c_str(), 100, 0, 3000}, "sumRms_cosine", "simu_energy_w_pathc");
        h_BGOrec_sumRms_bin_cosine2D[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                        .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                        .Histo2D<double, double>({(std::string("h_BGOrec_sumRms_bin_cosine2D_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); sumRms [mm]")).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
        h_BGOrec_sumRms_flast_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                     .Histo2D<double, double, double>({(std::string("h_BGOrec_sumRms_flast_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms vs F_{last} - bin ") + std::to_string(bin_idx) + std::string("; sumRms [mm]; F_{last}")).c_str(), (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
        h_BGOrec_sumRms_flast_13_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                        .Histo2D<double, double, double>({(std::string("h_BGOrec_sumRms_flast_13_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms vs F_{13} - bin ") + std::to_string(bin_idx) + std::string("sumRms [mm]; F_{13}")).c_str(), (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
        h_BGOrec_max_energy_ratio[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                     .Define("max_energy_ratio", GetMaxEnergyRatio, {"fracLayer"})
                                                     .Histo1D<double, double>({(std::string("h_BGOrec_max_energy_ratio_bin_") + std::to_string(bin_idx)).c_str(), (std::string("max energy ratio - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "max_energy_ratio", "simu_energy_w_pathc");
        h_BGOrec_ratio_last[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                               .Histo1D<double, double>({(std::string("h_BGOrec_ratio_last_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 0.01}, "fracLast", "simu_energy_w_pathc");
        h_BGOrec_ratio_13[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                             .Define("energy_ratio_13", "fracLayer[13]")
                                             .Histo1D<double, double>({(std::string("h_BGOrec_ratio_13_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 13th Layer - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, .1}, "energy_ratio_13", "simu_energy_w_pathc");
        h_BGOrec_ratio_last_cosine[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                      .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                      .Histo2D<double, double, double>({(std::string("h_BGOrec_ratio_last_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}")).c_str(), 100, 0, 1, 100, 0, 0.01}, "bgorec_cosine", "fracLast", "simu_energy_w_pathc");
        h_BGOrec_last_layer[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                               .Histo1D<int, double>({(std::string("h_BGOrec_last_layer_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Last BGO layer - bin ") + std::to_string(bin_idx)).c_str(), 14, 0, 13}, "lastBGOLayer", "simu_energy_w_pathc");
        h_BGOrec_hits[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                         .Histo1D<int, double>({(std::string("h_BGOrec_hits_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO hits - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1000}, "nBGOentries", "simu_energy_w_pathc");
        h_xtrl_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                      .Histo1D<double, double>({(std::string("h_xtrl_bin_") + std::to_string(bin_idx)).c_str(), (std::string("XTRL - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 150}, "xtrl", "simu_energy_w_pathc");

        _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
            .Foreach([&h_BGOrec_shower_profile, bin_idx](std::vector<double> &energy_frac, double energy_w) { 
                for (int lidx=0; lidx<DAMPE_bgo_nLayers; ++lidx)
                   h_BGOrec_shower_profile[bin_idx - 1]->Fill(lidx, energy_frac[lidx], energy_w); }, {"fracLayer", "simu_energy_w_pathc"});

        _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
            .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
            .Filter("bgorec_cosine<0.9")
            .Foreach([&h_BGOrec_shower_profile_upto_09, bin_idx](std::vector<double> &energy_frac, double energy_w) { 
                for (int lidx=0; lidx<DAMPE_bgo_nLayers; ++lidx)
                    h_BGOrec_shower_profile_upto_09[bin_idx - 1]->Fill(lidx, energy_frac[lidx], energy_w); }, {"fracLayer", "simu_energy_w_pathc"});

        _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
            .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
            .Filter("bgorec_cosine>0.9")
            .Foreach([&h_BGOrec_shower_profile_from_09, bin_idx](std::vector<double> &energy_frac, double energy_w) { 
                for (int lidx=0; lidx<DAMPE_bgo_nLayers; ++lidx)
                    h_BGOrec_shower_profile_from_09[bin_idx - 1]->Fill(lidx, energy_frac[lidx], energy_w); }, {"fracLayer", "simu_energy_w_pathc"});

        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto GetLayerComponent = [=](std::vector<double> &value_layer) -> double { return value_layer[ly]; };
            h_BGOrec_rms_layer[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                      .Define("rms_layer", GetLayerComponent, {"rmsLayer"})
                                                      .Histo1D<double, double>({(std::string("h_BGOrec_rms_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("RMS layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 500}, "rms_layer", "simu_energy_w_pathc");
            h_BGOrec_energy_ratio_layer[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                               .Define("energy_ratio_layer", GetLayerComponent, {"fracLayer"})
                                                               .Histo1D<double, double>({(std::string("h_BGOrec_energy_ratio_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_layer", "simu_energy_w_pathc");
            h_BGOrec_energy_ratio_1R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                            .Define("energy_1R_layer", GetLayerComponent, {"energy_1R_radius"})
                                                            .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                            .Define("energy_ratio_1R_layer", "energy_1R_layer/energy_layer")
                                                            .Histo1D<double, double>({(std::string("h_BGOrec_energy_ratio_1R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 1MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_1R_layer", "simu_energy_w_pathc");
            h_BGOrec_energy_ratio_2R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                            .Define("energy_2R_layer", GetLayerComponent, {"energy_2R_radius"})
                                                            .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                            .Define("energy_ratio_2R_layer", "energy_2R_layer/energy_layer")
                                                            .Histo1D<double, double>({(std::string("h_BGOrec_energy_ratio_2R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 2MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_2R_layer", "simu_energy_w_pathc");
            h_BGOrec_energy_ratio_3R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                            .Define("energy_3R_layer", GetLayerComponent, {"energy_3R_radius"})
                                                            .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                            .Define("energy_ratio_3R_layer", "energy_3R_layer/energy_layer")
                                                            .Histo1D<double, double>({(std::string("h_BGOrec_energy_ratio_3R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 3MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_3R_layer", "simu_energy_w_pathc");
            h_BGOrec_energy_ratio_5R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                            .Define("energy_5R_layer", GetLayerComponent, {"energy_5R_radius"})
                                                            .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                            .Define("energy_ratio_5R_layer", "energy_5R_layer/energy_layer")
                                                            .Histo1D<double, double>({(std::string("h_BGOrec_energy_ratio_5R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 5MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_5R_layer", "simu_energy_w_pathc");
        }
    }

    // Extract Simu histos
    auto h_simu_energy = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                             .Histo1D<double, double>({"h_simu_energy", "Simu energy", energy_nbins, &energy_binning[0]}, "simu_energy_gev", "simu_energy_w_pathc");
    auto h_simu_energy_w = _fr_bgo_analysis.Histo1D<double>({"h_simu_energy_w", "Simu energy weight", 100, 0, 1}, "simu_energy_w_pathc");
    auto h_energy_diff = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                             .Define("raw_energy_gev", "energy * 0.001")
                             .Define("energy_diff", "simu_energy_gev - raw_energy_gev")
                             .Histo1D<double, double>({"h_energy_diff", "Simu vs Raw Reco BGO energy: Simu Energy - Raw Energy (GeV); counts", 100, 0, 100}, "energy_diff", "simu_energy_w_pathc");
    auto h_energy_diff_corr = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                                  .Define("corr_energy_gev", "energy_corr * 0.001")
                                  .Define("energy_diff", "simu_energy_gev - corr_energy_gev")
                                  .Histo1D<double, double>({"h_energy_diff_corr", "Simu vs Corrected Reco BGO energy: Simu Energy - Corrected Energy (GeV); counts", 100, -100, 100}, "energy_diff", "simu_energy_w_pathc");
    auto h_energy_diff2D = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                               .Define("energy_ratio", "(energy - simu_energy)/simu_energy")
                               .Histo2D<double, double, double>({"h_energy_diff2D", "Energy Ratio; Real Energy (GeV); (Raw - Simu)/Simu", energy_nbins, &energy_binning[0], (int)energy_ratio_bins.size() - 1, &(energy_ratio_bins[0])}, "simu_energy_gev", "energy_ratio", "simu_energy_w_pathc");
    auto h_energy_diff2D_corr = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Define("energy_ratio", "(energy_corr-simu_energy)/simu_energy")
                                    .Histo2D<double, double, double>({"h_energy_diff2D_corr", "Energy Ratio; Real Energy (GeV); (Corr - Simu)/Simu", energy_nbins, &energy_binning[0], (int)energy_ratio_bins.size() - 1, &(energy_ratio_bins[0])}, "simu_energy_gev", "energy_ratio", "simu_energy_w_pathc");
    auto h_energy_unfold = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                               .Define("raw_energy_gev", "energy * 0.001")
                               .Histo2D<double, double, double>({"h_energy_unfold", "Energy Unfolding Matrix; Real Energy (GeV); Raw Energy (GeV)", energy_nbins, &energy_binning[0], energy_nbins, &energy_binning[0]}, "simu_energy_gev", "raw_energy_gev", "simu_energy_w_pathc");
    auto h_energy_unfold_corr = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Histo2D<double, double, double>({"h_energy_unfold_corr", "Energy Unfolding Matrix; Real Energy (GeV); Corr Energy (GeV)", energy_nbins, &energy_binning[0], energy_nbins, &energy_binning[0]}, "simu_energy_gev", "corr_energy_gev", "simu_energy_w_pathc");
    auto h_simu_slopeX = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_slopeX", "Simu Slope X", 100, -10, 10}, "simuSlopeX", "simu_energy_w_pathc");
    auto h_simu_slopeY = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_slopeY", "Simu Slope Y", 100, -10, 10}, "simuSlopeY", "simu_energy_w_pathc");
    auto h_simu_interceptX = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_interceptX", "Simu Intercept X", 100, -500, 500}, "simuInterceptX", "simu_energy_w_pathc");
    auto h_simu_interceptY = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_interceptY", "Simu Intercept Y", 100, -500, 500}, "simuInterceptY", "simu_energy_w_pathc");

    // Extract STK histos
    auto h_stk_cosine = _fr_stk_analysis.Histo1D({"h_stk_cosine", "h_stk_cosine", 100, 0, 1}, "STK_bestTrack_costheta");
    auto h_stk_slopeX = _fr_stk_analysis.Histo1D<double, double>({"h_stk_slopeX", "Simu Slope X", 200, -10, 10}, "STK_bestTrack_slopeX", "simu_energy_w_pathc");
    auto h_stk_slopeY = _fr_stk_analysis.Histo1D<double, double>({"h_stk_slopeY", "Simu Slope Y", 200, -10, 10}, "STK_bestTrack_slopeY", "simu_energy_w_pathc");
    auto h_stk_interceptX = _fr_stk_analysis.Histo1D<double, double>({"h_stk_interceptX", "Simu Intercept X", 100, -500, 500}, "STK_bestTrack_interceptX", "simu_energy_w_pathc");
    auto h_stk_interceptY = _fr_stk_analysis.Histo1D<double, double>({"h_stk_interceptY", "Simu Intercept Y", 100, -500, 500}, "STK_bestTrack_interceptY", "simu_energy_w_pathc");

    // Extract PSD charge histos
    auto h_psd_chargeX = _fr_psd_charge_analysis.Histo1D<double, double > ({"h_psd_chargeX", "PSD Charge X", 100, 0, 20}, "PSD_chargeX", "simu_energy_w_pathc");
    auto h_psd_chargeY = _fr_psd_charge_analysis.Histo1D<double, double > ({"h_psd_chargeY", "PSD Charge Y", 100, 0, 20}, "PSD_chargeY", "simu_energy_w_pathc");
    auto h_psd_charge = _fr_psd_charge_analysis.Define("psd_charge", [](double psd_chargeX, double psd_chargeY) { return 0.5 * (psd_chargeX + psd_chargeY); }, {"PSD_chargeX", "PSD_chargeY"})
                            .Histo1D<double, double>({"h_psd_charge", "PSD Charge", 100, 0, 20}, "psd_charge", "simu_energy_w_pathc");
    auto h_psd_charge2D = _fr_psd_charge_analysis.Histo2D<double, double, double>({"h_psd_charge2D", "PSD Charge", 100, 0, 20, 100, 0, 20}, "PSD_chargeX", "PSD_chargeY", "simu_energy_w_pathc");

    auto h_psd_selected_chargeX = _fr_psd_charge_analysis.Filter("evtfilter_psd_charge_cut==true")
                                      .Histo1D<double, double>({"h_psd_selected_chargeX", "PSD Charge X", 100, 0, 20}, "PSD_chargeX", "simu_energy_w_pathc");
    auto h_psd_selected_chargeY = _fr_psd_charge_analysis.Filter("evtfilter_psd_charge_cut==true")
                                      .Histo1D<double, double > ({"h_psd_selected_chargeY", "PSD Charge Y", 100, 0, 20}, "PSD_chargeY", "simu_energy_w_pathc");
    auto h_psd_selected_charge = _fr_psd_charge_analysis.Filter("evtfilter_psd_charge_cut==true")
                                     .Define("psd_charge", [](double psd_chargeX, double psd_chargeY) { return 0.5 * (psd_chargeX + psd_chargeY); }, {"PSD_chargeX", "PSD_chargeY"})
                                     .Histo1D<double, double>({"h_psd_selected_charge", "PSD Charge", 100, 0, 20}, "psd_charge", "simu_energy_w_pathc");
    auto h_psd_selected_charge2D = _fr_psd_charge_analysis.Filter("evtfilter_psd_charge_cut==true")
                                       .Histo2D<double, double, double>({"h_psd_selected_charge2D", "PSD Charge", 100, 0, 20, 100, 0, 20}, "PSD_chargeX", "PSD_chargeY", "simu_energy_w_pathc");

    // Extract STK charge histos
    auto h_stk_chargeX = _fr_stk_charge_analysis.Histo1D<double, double > ({"h_stk_chargeX", "STK Charge X", 100, 0, 20}, "STK_chargeX", "simu_energy_w_pathc");
    auto h_stk_chargeY = _fr_stk_charge_analysis.Histo1D<double, double > ({"h_stk_chargeY", "STK Charge Y", 100, 0, 20}, "STK_chargeY", "simu_energy_w_pathc");
    auto h_stk_charge = _fr_stk_charge_analysis.Define("stk_charge", [](double stk_chargeX, double stk_chargeY) { return 0.5 * (stk_chargeX + stk_chargeY); }, {"STK_chargeX", "STK_chargeY"})
                            .Histo1D<double, double>({"h_stk_charge", "STK Charge", 100, 0, 20}, "stk_charge", "simu_energy_w_pathc");
    auto h_stk_charge2D = _fr_stk_charge_analysis.Histo2D<double, double, double>({"h_stk_charge2D", "STK Charge", 100, 0, 20, 100, 0, 20}, "STK_chargeX", "STK_chargeY", "simu_energy_w_pathc");

    auto h_stk_selected_chargeX = _fr_stk_charge_analysis.Filter("evtfilter_stk_charge_cut==true")
                                      .Histo1D<double, double>({"h_stk_selected_chargeX", "STK Charge X", 100, 0, 20}, "STK_chargeX", "simu_energy_w_pathc");
    auto h_stk_selected_chargeY = _fr_stk_charge_analysis.Filter("evtfilter_stk_charge_cut==true")
                                      .Histo1D<double, double > ({"h_stk_selected_chargeY", "STK Charge Y", 100, 0, 20}, "STK_chargeY", "simu_energy_w_pathc");
    auto h_stk_selected_charge = _fr_stk_charge_analysis.Filter("evtfilter_stk_charge_cut==true")
                                     .Define("stk_charge", [](double stk_chargeX, double stk_chargeY) { return 0.5 * (stk_chargeX + stk_chargeY); }, {"STK_chargeX", "STK_chargeY"})
                                     .Histo1D<double, double>({"h_stk_selected_charge", "STK Charge", 100, 0, 20}, "stk_charge", "simu_energy_w_pathc");
    auto h_stk_selected_charge2D = _fr_stk_charge_analysis.Filter("evtfilter_stk_charge_cut==true")
                                       .Histo2D<double, double, double>({"h_stk_selected_charge2D", "STK Charge", 100, 0, 20, 100, 0, 20}, "STK_chargeX", "STK_chargeY", "simu_energy_w_pathc");

    // Extract NUD hisos
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_NUD_adc (DAMPE_NUD_channels); 

    for (int channel = 0; channel < DAMPE_NUD_channels; ++channel)
        h_NUD_adc[channel] = _fr_preselected.Define("nud_adc_channel", [channel](std::vector<double> nud_adc) {return nud_adc[channel];}, {"NUD_ADC"})
                                            .Histo1D<double, double>({(std::string("h_NUD_adc_") + std::to_string(channel)).c_str(), (std::string("NUD ADC - channel ") + std::to_string(channel)).c_str(), 100, 0, 1000}, "nud_adc_channel", "simu_energy_w_pathc");
    /*
    auto h_NUD_total_adc = _fr_preselected.Histo1D<double, double>({"h_NUD_total_adc", "NUD Total ADC", 100, 0, 10000}, "NUD_total_ADC", "simu_energy_w_pathc");
    auto h_NUD_max_adc = _fr_preselected.Histo1D<double, double>({"h_NUD_max_adc", "NUD Max ADC", 100, 0, 1000}, "NUD_max_ADC", "simu_energy_w_pathc");
    auto h_NUD_max_channel = _fr_preselected.Histo1D<double, double>({"h_NUD_max_channel", "NUD Max Channel", 3, 0, 3}, "NUD_max_channel_ID", "simu_energy_w_pathc");
    */

    if (_VERBOSE)
        std::cout << "Writing to disk... [" << outputPath << "]" << std::endl;

    TFile *outfile = TFile::Open(outputPath.c_str(), "RECREATE");

    auto simu_dir = outfile->mkdir("Simu");
    simu_dir->cd();

    h_simu_energy->Write();
    h_simu_energy_w->Write();
    h_energy_diff->Write();
    h_energy_diff_corr->Write();
    h_energy_diff2D->Write();
    h_energy_diff2D_corr->Write();
    h_energy_unfold->Write();
    h_energy_unfold_corr->Write();
    h_simu_slopeX->Write();
    h_simu_slopeY->Write();
    h_simu_interceptX->Write();
    h_simu_interceptY->Write();

    auto stk_dir = outfile->mkdir("STK");
    stk_dir->cd();

    h_stk_cosine->Write();
    h_stk_slopeX->Write();
    h_stk_slopeY->Write();
    h_stk_interceptX->Write();
    h_stk_interceptY->Write();

    h_stk_chargeX->Write();
    h_stk_chargeY->Write();
    h_stk_charge->Write();
    h_stk_charge2D->Write();
    h_stk_selected_chargeX->Write();
    h_stk_selected_chargeY->Write();
    h_stk_selected_charge->Write();
    h_stk_selected_charge2D->Write();

    auto psd_dir = outfile->mkdir("PSD");
    psd_dir->cd();

    h_psd_chargeX->Write();
    h_psd_chargeY->Write();
    h_psd_charge->Write();
    h_psd_charge2D->Write();
    h_psd_selected_chargeX->Write();
    h_psd_selected_chargeY->Write();
    h_psd_selected_charge->Write();
    h_psd_selected_charge2D->Write();

    auto bgo_dir = outfile->mkdir("BGO");
    bgo_dir->cd();

    h_BGOrec_raw_energy->Write();
    h_BGOrec_corr_energy->Write();
    h_BGOrec_cosine->Write();
    h_BGOrec_cosine_preseltion->Write();
    h_BGOrec_slopeX->Write();
    h_BGOrec_slopeY->Write();
    h_BGOrec_interceptX->Write();
    h_BGOrec_interceptY->Write();
    h_BGOrec_sumRms_bin_cosine2D_20_100->Write();
    h_BGOrec_sumRms_bin_cosine2D_100_250->Write();
    h_BGOrec_sumRms_bin_cosine2D_250_500->Write();
    h_BGOrec_sumRms_bin_cosine2D_500_1000->Write();
    h_BGOrec_sumRms_bin_cosine2D_1000_3000->Write();
    h_BGOrec_sumRms_bin_cosine2D_3000_10000->Write();
    h_BGOrec_sumRms_bin_cosine2D_10000_20000->Write();
    h_BGOrec_sumRms_flast_20_100->Write();
    h_BGOrec_sumRms_flast_100_250->Write();
    h_BGOrec_sumRms_flast_250_500->Write();
    h_BGOrec_sumRms_flast_500_1000->Write();
    h_BGOrec_sumRms_flast_1000_3000->Write();
    h_BGOrec_sumRms_flast_3000_10000->Write();
    h_BGOrec_sumRms_flast_10000_20000->Write();
    h_BGOrec_sumRms_flast_13_20_100->Write();
    h_BGOrec_sumRms_flast_13_100_250->Write();
    h_BGOrec_sumRms_flast_13_250_500->Write();
    h_BGOrec_sumRms_flast_13_500_1000->Write();
    h_BGOrec_sumRms_flast_13_1000_3000->Write();
    h_BGOrec_sumRms_flast_13_3000_10000->Write();
    h_BGOrec_sumRms_flast_13_10000_20000->Write();
    h_xtrl_energy_int->Write();
    h_xtrl->Write();

    for (int bidx = 0; bidx < energy_nbins; ++bidx)
    {
        auto tmp_dir_name = std::string("BGO/energybin_") + std::to_string(bidx + 1);
        outfile->mkdir(tmp_dir_name.c_str());
        outfile->cd(tmp_dir_name.c_str());

        h_BGOrec_cosine_bin[bidx]->Write();
        h_BGOrec_sumRms_bin[bidx]->Write();
        h_BGOrec_sumRms_mean_bin[bidx]->Write();
        h_BGOrec_sumRms_bin_weight[bidx]->Write();
        h_BGOrec_sumRms_sumRms_weight[bidx]->Write();
        h_BGOrec_sumRms_bin_cosine[bidx]->Write();
        h_BGOrec_sumRms_bin_cosine2D[bidx]->Write();
        h_BGOrec_sumRms_flast_bin[bidx]->Write();
        h_BGOrec_sumRms_flast_13_bin[bidx]->Write();
        h_BGOrec_max_energy_ratio[bidx]->Write();
        h_BGOrec_ratio_last[bidx]->Write();
        h_BGOrec_ratio_13[bidx]->Write();
        h_BGOrec_ratio_last_cosine[bidx]->Write();
        h_BGOrec_last_layer[bidx]->Write();
        h_BGOrec_hits[bidx]->Write();
        h_xtrl_bin[bidx]->Write();
        h_BGOrec_shower_profile[bidx]->Write();
        h_BGOrec_shower_profile_upto_09[bidx]->Write();
        h_BGOrec_shower_profile_from_09[bidx]->Write();

        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            h_BGOrec_rms_layer[bidx][ly]->Write();
            h_BGOrec_energy_ratio_layer[bidx][ly]->Write();
            h_BGOrec_energy_ratio_1R[bidx][ly]->Write();
            h_BGOrec_energy_ratio_2R[bidx][ly]->Write();
            h_BGOrec_energy_ratio_3R[bidx][ly]->Write();
            h_BGOrec_energy_ratio_5R[bidx][ly]->Write();
        }
    }

    auto nud_dir = outfile->mkdir("NUD");
    nud_dir->cd();

    for (auto& _elm : h_NUD_adc)
        _elm->Write();
    /*
    h_NUD_total_adc->Write();
    h_NUD_max_adc->Write();
    h_NUD_max_channel->Write();
    */
   
    outfile->Close();
}