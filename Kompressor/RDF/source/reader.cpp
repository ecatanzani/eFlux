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
        std::cout << "Total number of events: " << _entries << "\n\n";
        std::cout << "Reading events..." << std::endl;
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
    ROOT::EnableImplicitMT(nthreads);
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
    auto good_for_bgo_analysis = [](bool out_energy, bool trigger, bool bgo_reco) { return (!out_energy) * trigger * bgo_reco; };
    auto preselection = [](bool all_cut) { return all_cut; };
    // Create the RDF for the BGO analysis
    auto _fr_bgo_analysis = _fr_weight_patch.Filter(good_for_bgo_analysis, {"evtfilter_out_energy_range", "evtfilter_evt_triggered", "evtfilter_correct_bgo_reco"});
    //Create the RDF for the preselected events
    auto _fr_preselected = _fr_bgo_analysis.Filter(preselection, {"evtfilter_all_cut"});

    // Extract BGO histos
    auto sumRms_bins = createLogBinning(10, 2e+3, 1e+2);
    auto xtrl_bins = createLinearBinning(0, 150, 1e+2);
    auto cosine_bins = createLinearBinning(0, 1, 1e+2);
    auto flast_binning = createLogBinning(1e-5, 2e-1, 1e+2);

    auto h_BGOrec_raw_energy = _fr_bgo_analysis.Define("raw_energy_gev", "energy * 0.001")
                                   .Histo1D<double>({"h_BGOrec_raw_energy", "BGO Raw energy", energy_nbins, &energy_binning[0]}, "raw_energy_gev", "simu_energy_w_pathc");
    auto h_BGOrec_corr_energy = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Histo1D<double>({"h_BGOrec_corr_energy", "BGO Corr energy", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "simu_energy_w_pathc");
    auto h_BGOrec_cosine = _fr_bgo_analysis.Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                               .Histo1D<double>({"h_BGOrec_cosine", "BGO Reco cosine", 100, 0, 1}, "bgorec_cosine", "simu_energy_w_pathc");
    auto h_BGOrec_cosine_preseltion = _fr_preselected.Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                          .Histo1D<double>({"h_BGOrec_cosine_preseltion", "BGO Reco cosine - Preselection", 100, 0, 1}, "bgorec_cosine", "simu_energy_w_pathc");

    auto h_sumRms_bin_cosine2D_20_100 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                            .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                            .Histo2D<double>({"h_sumRms_bin_cosine2D_20_100", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_sumRms_bin_cosine2D_100_250 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                             .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                             .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                             .Histo2D<double>({"h_sumRms_bin_cosine2D_100_250", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_sumRms_bin_cosine2D_250_500 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                             .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                             .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                             .Histo2D<double>({"h_sumRms_bin_cosine2D_250_500", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_sumRms_bin_cosine2D_500_1000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                              .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                              .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                              .Histo2D<double>({"h_sumRms_bin_cosine2D_500_1000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_sumRms_bin_cosine2D_1000_3000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                               .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                               .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                               .Histo2D<double>({"h_sumRms_bin_cosine2D_1000_3000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_sumRms_bin_cosine2D_3000_10000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                .Histo2D<double>({"h_sumRms_bin_cosine2D_3000_10000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    auto h_sumRms_bin_cosine2D_10000_20000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                 .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                                 .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                 .Histo2D<double>({"h_sumRms_bin_cosine2D_10000_10000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
    
    auto h_sumRms_flast_20_100 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                     .Histo2D<double>({"h_sumRms_flast_20_100", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_sumRms_flast_100_250 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                     .Histo2D<double>({"h_sumRms_flast_100_250", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_sumRms_flast_250_500 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                     .Histo2D<double>({"h_sumRms_flast_250_500", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_sumRms_flast_500_1000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                     .Histo2D<double>({"h_sumRms_flast_500_1000", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_sumRms_flast_1000_3000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                     .Histo2D<double>({"h_sumRms_flast_1000_3000", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_sumRms_flast_3000_10000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                     .Histo2D<double>({"h_sumRms_flast_3000_10000", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
    auto h_sumRms_flast_10000_20000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                     .Histo2D<double>({"h_sumRms_flast_10000_20000", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");

    auto h_sumRms_flast_13_20_100 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                     .Histo2D<double>({"h_sumRms_flast_13_20_100", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_sumRms_flast_13_100_250 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                     .Histo2D<double>({"h_sumRms_flast_13_100_250", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_sumRms_flast_13_250_500 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                     .Histo2D<double>({"h_sumRms_flast_13_250_500", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_sumRms_flast_13_500_1000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                     .Histo2D<double>({"h_sumRms_flast_13_500_1000", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_sumRms_flast_13_1000_3000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                     .Histo2D<double>({"h_sumRms_flast_13_1000_3000", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_sumRms_flast_13_3000_10000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                     .Histo2D<double>({"h_sumRms_flast_13_3000_10000", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
    auto h_sumRms_flast_13_10000_20000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                     .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                     .Histo2D<double>({"h_sumRms_flast_13_10000_20000", "sumRms vs flast - bin ", (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");

    auto h_xtrl_energy_int = _fr_bgo_analysis.Histo1D<double>({"h_xtrl_energy_int", "XTRL", 100, 0, 150}, "xtrl", "simu_energy_w_pathc");
    auto h_xtrl = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                      .Histo2D<double>({"h_xtrl", "XTRL", energy_nbins, &energy_binning[0], (int)xtrl_bins.size() - 1, &xtrl_bins[0]}, "corr_energy_gev", "xtrl", "simu_energy_w_pathc");

    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_cosine_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_sumRms_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_sumRms_bin_weight(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_sumRms_bin_cosine(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_sumRms_bin_cosine2D(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_sumRms_flast_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_sumRms_flast_13_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_max_energy_ratio(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_ratio_last(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_ratio_13(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_ratio_last_cosine(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_last_layer(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_hits(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_shower_profile(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_shower_profile_upto_09(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_shower_profile_from_09(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_xtrl_bin(energy_nbins);
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_rms_layer(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_energy_ratio_layer(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_energy_ratio_1R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_energy_ratio_2R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_energy_ratio_3R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_energy_ratio_5R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));

    auto computeSumRmsWeight = [=](std::vector<double> &energy_layer, std::vector<double> &rms_layer, double raw_energy) {
        double sumRms_weight = 0;
        for (int ly = 0; ly<DAMPE_bgo_nLayers; ++ly)
            sumRms_weight += energy_layer[ly]*rms_layer[ly];
        return sumRms_weight/raw_energy; };
    auto GetMaxEnergyRatio = [](std::vector<double> &energy_fraction) -> double { 
        auto it_max = std::max_element(energy_fraction.begin(), energy_fraction.end());
        return energy_fraction[std::distance(energy_fraction.begin(), it_max)]; };

    std::vector<int> idx_layer (DAMPE_bgo_nLayers);
    std::iota (std::begin(idx_layer), std::end(idx_layer), 0);

    // Loop over energy bins
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        auto bin_filter = [=](int energy_bin) -> bool { return energy_bin == bin_idx; };
        h_BGOrec_cosine_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                               .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                               .Histo1D<double>({(std::string("h_BGOrec_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Reco cosine - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "bgorec_cosine", "simu_energy_w_pathc");
        h_sumRms_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                        .Histo1D<double>({(std::string("h_sumRms_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 3000}, "sumRms", "simu_energy_w_pathc");
        h_sumRms_bin_weight[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                               .Define("sumRms_w", computeSumRmsWeight, {"eLayer", "rmsLayer", "energy"})
                                               .Histo1D<double>({(std::string("h_sumRms_weight_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms weighted - bin ") + std::to_string(bin_idx)).c_str(), 50, 0, 100}, "sumRms_w", "simu_energy_w_pathc");
        h_sumRms_bin_cosine[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                               .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                               .Define("sumRms_cosine", "sumRms/bgorec_cosine")
                                               .Histo1D<double>({(std::string("h_sumRms_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 3000}, "sumRms_cosine", "simu_energy_w_pathc");
        h_sumRms_bin_cosine2D[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                 .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                 .Histo2D<double>({(std::string("h_sumRms_bin_cosine2D_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx)).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_bins.size() - 1, &sumRms_bins[0]}, "bgorec_cosine", "sumRms", "simu_energy_w_pathc");
        h_sumRms_flast_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                              .Histo2D<double>({(std::string("h_sumRms_flast_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms vs flast - bin ") + std::to_string(bin_idx)).c_str(), (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast", "simu_energy_w_pathc");
        h_sumRms_flast_13_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                 .Histo2D<double>({(std::string("h_sumRms_flast_13_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms vs flast 13th layer - bin ") + std::to_string(bin_idx)).c_str(), (int)sumRms_bins.size() - 1, &sumRms_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "sumRms", "fracLast_13", "simu_energy_w_pathc");
        h_max_energy_ratio[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                              .Define("max_energy_ratio", GetMaxEnergyRatio, {"fracLayer"})
                                              .Histo1D<double>({(std::string("h_max_energy_ratio_bin_") + std::to_string(bin_idx)).c_str(), (std::string("max energy ratio - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "max_energy_ratio", "simu_energy_w_pathc");
        h_ratio_last[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                        .Histo1D<double>({(std::string("h_ratio_last_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 0.01}, "fracLast", "simu_energy_w_pathc");
        h_ratio_13[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                      .Define("energy_ratio_13", "fracLayer[13]")
                                      .Histo1D<double>({(std::string("h_ratio_13_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 13th Layer - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, .1}, "energy_ratio_13", "simu_energy_w_pathc");
        h_ratio_last_cosine[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                               .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                               .Histo2D<double>({(std::string("h_ratio_last_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer vs cosine - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1, 100, 0, 0.01}, "bgorec_cosine", "fracLast", "simu_energy_w_pathc");
        h_last_layer[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                        .Histo1D<double>({(std::string("h_last_layer_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Last BGO layer - bin ") + std::to_string(bin_idx)).c_str(), 14, 0, 13}, "lastBGOLayer", "simu_energy_w_pathc");
        h_hits[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                  .Histo1D<double>({(std::string("h_hits_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO hits - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 4000}, "nBGOentries", "simu_energy_w_pathc");
        h_xtrl_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                      .Histo1D<double>({(std::string("h_xtrl_bin_") + std::to_string(bin_idx)).c_str(), (std::string("XTRL - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 150}, "xtrl", "simu_energy_w_pathc");
        
        auto get_frac_comp = [=] (std::vector<double>& vec) -> double {for (auto& _elm : vec) return _elm;};
        
        /*
        h_shower_profile[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                .Define("layer_index", [&idx_layer] { for (auto& _elm : idx_layer) return _elm; })
                                                .Define("energy_ratio_layer", get_frac_comp, {"fracLayer"})
                                                .Histo2D<int, double, double>({(std::string("h_shower_profile_bin") + std::to_string(bin_idx)).c_str(), (std::string("BGO Shower Profile - bin ") + std::to_string(bin_idx)).c_str(), 13, 0, 13, 100, 0, 0.01}, "layer_index", "energy_ratio_layer", "simu_energy_w_pathc");
        h_shower_profile_upto_09[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                    .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                    .Filter("bgorec_cosine <= 0.9")
                                                    .Define("layer_index", [=] { return ly; })
                                                    .Define("energy_ratio_layer", GetLayerComponent, {"fracLayer"})
                                                    .Histo2D<double>({(std::string("h_shower_profile_upto_09_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Shower Profile - bin ") + std::to_string(bin_idx)).c_str(), 14, 0, 13, 100, 0, 0.01}, "layer_index", "energy_ratio_layer", "simu_energy_w_pathc");
        h_shower_profile_from_09[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                    .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                    .Filter("bgorec_cosine > 0.9")
                                                    .Define("layer_index", [=] { return ly; })
                                                    .Define("energy_ratio_layer", GetLayerComponent, {"fracLayer"})
                                                    .Histo2D<double>({(std::string("h_shower_profile_from_09_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Shower Profile - bin ") + std::to_string(bin_idx)).c_str(), 14, 0, 13, 100, 0, 0.01}, "layer_index", "energy_ratio_layer", "simu_energy_w_pathc");
        */

        // Loop over BGO layers
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto GetLayerComponent = [=](std::vector<double> &value_layer) -> double { return value_layer[ly]; };
            h_rms_layer[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                               .Define("rms_layer", GetLayerComponent, {"rmsLayer"})
                                               .Histo1D<double>({(std::string("h_rms_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("RMS layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 500}, "rms_layer", "simu_energy_w_pathc");
            h_energy_ratio_layer[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                        .Define("energy_ratio_layer", GetLayerComponent, {"fracLayer"})
                                                        .Histo1D<double>({(std::string("h_energy_ratio_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_layer", "simu_energy_w_pathc");
            h_energy_ratio_1R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                     .Define("energy_1R_layer", GetLayerComponent, {"energy_1R_radius"})
                                                     .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                     .Define("energy_ratio_1R_layer", "energy_1R_layer/energy_layer")
                                                     .Histo1D<double>({(std::string("h_energy_ratio_1R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 1MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_1R_layer", "simu_energy_w_pathc");
            h_energy_ratio_2R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                     .Define("energy_2R_layer", GetLayerComponent, {"energy_2R_radius"})
                                                     .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                     .Define("energy_ratio_2R_layer", "energy_2R_layer/energy_layer")
                                                     .Histo1D<double>({(std::string("h_energy_ratio_2R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 2MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_2R_layer", "simu_energy_w_pathc");
            h_energy_ratio_3R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                     .Define("energy_3R_layer", GetLayerComponent, {"energy_3R_radius"})
                                                     .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                     .Define("energy_ratio_3R_layer", "energy_3R_layer/energy_layer")
                                                     .Histo1D<double>({(std::string("h_energy_ratio_3R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 3MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_3R_layer", "simu_energy_w_pathc");
            h_energy_ratio_5R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                     .Define("energy_5R_layer", GetLayerComponent, {"energy_5R_radius"})
                                                     .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                     .Define("energy_ratio_5R_layer", "energy_5R_layer/energy_layer")
                                                     .Histo1D<double>({(std::string("h_energy_ratio_5R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 5MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_5R_layer", "simu_energy_w_pathc");
        }
    }

    //auto h_stk_cosine = _data_fr.Filter("STK_bestTrack_costheta > 0").Histo1D({"h_stk_cosine", "h_stk_cosine", 100, 0, 1}, "STK_bestTrack_costheta");

    TFile *outfile = TFile::Open(outputPath.c_str(), "RECREATE");
    //h_stk_cosine->Write();

    auto bgo_dir = outfile->mkdir("BGO");
    bgo_dir->cd();
    h_BGOrec_raw_energy->Write();
    h_BGOrec_corr_energy->Write();
    h_BGOrec_cosine->Write();
    h_BGOrec_cosine_preseltion->Write();
    h_sumRms_bin_cosine2D_20_100->Write();
    h_sumRms_bin_cosine2D_100_250->Write();
    h_sumRms_bin_cosine2D_250_500->Write();
    h_sumRms_bin_cosine2D_500_1000->Write();
    h_sumRms_bin_cosine2D_1000_3000->Write();
    h_sumRms_bin_cosine2D_3000_10000->Write();
    h_sumRms_bin_cosine2D_10000_20000->Write();
    h_sumRms_flast_20_100->Write();
    h_sumRms_flast_100_250->Write();
    h_sumRms_flast_250_500->Write();
    h_sumRms_flast_500_1000->Write();
    h_sumRms_flast_1000_3000->Write();
    h_sumRms_flast_3000_10000->Write();
    h_sumRms_flast_10000_20000->Write();
    h_sumRms_flast_13_20_100->Write();
    h_sumRms_flast_13_100_250->Write();
    h_sumRms_flast_13_250_500->Write();
    h_sumRms_flast_13_500_1000->Write();
    h_sumRms_flast_13_1000_3000->Write();
    h_sumRms_flast_13_3000_10000->Write();
    h_sumRms_flast_13_10000_20000->Write();
    h_xtrl_energy_int->Write();
    h_xtrl->Write();

    for (int bidx = 0; bidx < energy_nbins; ++bidx)
    {
        auto tmp_dir_name = std::string("BGO/energybin_") + std::to_string(bidx + 1);
        outfile->mkdir(tmp_dir_name.c_str());
        outfile->cd(tmp_dir_name.c_str());

        h_BGOrec_cosine_bin[bidx]->Write();
        h_sumRms_bin[bidx]->Write();
        h_sumRms_bin_weight[bidx]->Write();
        h_sumRms_bin_cosine[bidx]->Write();
        h_sumRms_bin_cosine2D[bidx]->Write();
        h_sumRms_flast_bin[bidx]->Write();
        h_sumRms_flast_13_bin[bidx]->Write();
        h_max_energy_ratio[bidx]->Write();
        h_ratio_last[bidx]->Write();
        h_ratio_13[bidx]->Write();
        h_ratio_last_cosine[bidx]->Write();
        h_last_layer[bidx]->Write();
        h_hits[bidx]->Write();
        /*
        h_shower_profile[bidx]->Write();
        h_shower_profile_upto_09[bidx]->Write();
        h_shower_profile_from_09[bidx]->Write();
        */
        h_xtrl_bin[bidx]->Write();

        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            h_rms_layer[bidx][ly]->Write();
            h_energy_ratio_layer[bidx][ly]->Write();
            h_energy_ratio_1R[bidx][ly]->Write();
            h_energy_ratio_2R[bidx][ly]->Write();
            h_energy_ratio_3R[bidx][ly]->Write();
            h_energy_ratio_5R[bidx][ly]->Write();
        }
    }

    outfile->Close();
}