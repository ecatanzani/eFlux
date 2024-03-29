#include "reader.h"
#include "binning.h"
#include "kompress.h"
#include "regularize.h"
#include "DAMPE_geo_structure.h"

#include "TF1.h"
#include "TFile.h"
#include "TVector3.h"
#include <ROOT/RDataFrame.hxx>

inline std::vector<std::shared_ptr<TH2D>> extractHistoFrorRegBH(const int bins, const char* path="/mnt/eFlux/Kompressor/protons_kompressor_0123_sumrms_reg_by_hand_integral_6.root") {
    TFile* infile = TFile::Open(path, "READ");
    if (infile->IsZombie()) {
        std::cerr << "\n\nError reading Reg Histo BH: [" << path << "]\n\n";
        exit(100);
    }

    std::vector<std::shared_ptr<TH2D>> histos;

    for (int idx=0; idx<bins; ++idx) {
        std::string histo_name = "Preselection/BGO/energybin_" + std::to_string(idx+1) + std::string("/h_BGOrec_ps_sumRms_cosine2D_bin_") + std::to_string(idx+1);
        auto tmp_histo = std::shared_ptr<TH2D>(static_cast<TH2D*>(infile->Get(histo_name.c_str())));
        tmp_histo->SetDirectory(0);
        histos.push_back(tmp_histo);
    }

    infile->Close();

    return histos;
}

void kompress(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    const std::string fit_tree_path,
    const double _entries,
    const std::string outputPath,
    const bool verbose,
    const unsigned int threads,
    const bool _mc) {

    // Extract the energy info
    auto energy_binning = _energy_config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;
    auto min_evt_energy = _energy_config->GetMinEvtEnergy();
    auto max_evt_energy = _energy_config->GetMaxEvtEnergy();
    // Plotting variables
    double _gev = 0.001;
    std::string sumRms_leaf = "sumRms";
    std::string fracLast_leaf = "fracLast";
    // Set the regularization options
    bool regularize_vars;
    fit_tree_path.empty() ? regularize_vars = false : regularize_vars = true;
    
    // Create the RDF
    ROOT::EnableImplicitMT(threads);
    ROOT::RDataFrame _data_fr(*evtch);
    
    auto GetEnergyBin = [=](double energy) -> int { 
        int bin_idx=0;
        for (; bin_idx<energy_nbins-1; ++bin_idx)
            if (energy * _gev >= energy_binning[bin_idx] && energy * _gev < energy_binning[bin_idx+1])
                break;
        return bin_idx+1; 
    };
    
    auto energyFilter = [&min_evt_energy, &max_evt_energy, &_gev](const double energy) -> bool {
        return energy*_gev >= min_evt_energy && energy*_gev <= max_evt_energy ? true : false;
    };

    auto simuEnergyWeight = [&energy_binning](const double simu_energy_w) -> double {
        return simu_energy_w * pow(energy_binning[0], 2);
    };

    // Create the RDF with energy bin and correct energy weight
    auto _fr_energy_filter = _data_fr.Filter(energyFilter, {"energy_corr"})
                                        .Define("energy_bin", GetEnergyBin, {"energy_corr"});
    
    _fr_energy_filter = _mc ? _fr_energy_filter.Define("simu_energy_w_corr", simuEnergyWeight, {"simu_energy_w"}) : _fr_energy_filter.Define("simu_energy_w_corr", "1.");

    // Regularize RDF
    std::vector<TF1> sumrms_fitfunc(energy_nbins);
    std::vector<TF1> sumrms_fitfunc_err(energy_nbins);
    std::vector<TF1> flast_fitfunc(energy_nbins);
    std::vector<TF1> flast_fitfunc_err(energy_nbins);

    auto regularize_sumrms = [&sumrms_fitfunc, &sumrms_fitfunc_err](double sumrms, int energy_bin, TVector3 bgodir) -> double {
        // Initialize regularized sumrms variable
        double reg_sumrms = sumrms;
        // Initialize BGO cosine from directrion
        double bgocosine = bgodir.CosTheta();
        // Regularize sumrms
        reg_sumrms -= sumrms_fitfunc[energy_bin - 1].Eval(bgocosine);
        reg_sumrms /= sumrms_fitfunc_err[energy_bin - 1].Eval(bgocosine);
        return reg_sumrms;
    };

    auto regularize_flast = [&flast_fitfunc, &flast_fitfunc_err](double flast, int energy_bin, TVector3 bgodir) -> double {
        // Initialize regularized sumrms variable
        double reg_flast = flast;
        // Initialize BGO cosine from directrion
        double bgocosine = bgodir.CosTheta();
        // Regularize sumrms
        reg_flast -= flast_fitfunc[energy_bin - 1].Eval(bgocosine);
        reg_flast /= flast_fitfunc_err[energy_bin - 1].Eval(bgocosine);
        return reg_flast;
    };

    if (regularize_vars) {
        if (verbose) std::cout << "\n\nINFO: Normalizing data facility has been activated\n";
        // Reading fit summary TTree
        load_tf1s(
            fit_tree_path,
            sumrms_fitfunc,
            sumrms_fitfunc_err,
            flast_fitfunc,
            flast_fitfunc_err,
            verbose);

        _fr_energy_filter = _fr_energy_filter.Define("sumRms_reg", regularize_sumrms, {"sumRms", "energy_bin", "BGOrec_trajectoryDirection2D"})
                                                .Define("fracLast_reg", regularize_flast, {"fracLast", "energy_bin", "BGOrec_trajectoryDirection2D"});

        sumRms_leaf = "sumRms_reg";
        fracLast_leaf = "fracLast_reg";
    }

    // Create the RDF for BGO analysis
    auto _fr_bgo_analysis = _fr_energy_filter.Filter("evtfilter_good_event==true");
    // Create the RDF for the preselected events
    auto _fr_preselected = _fr_bgo_analysis.Filter("evtfilter_all_cut==true");
    // Create the RDF for STK analysis
    auto _fr_stk_analysis = _fr_bgo_analysis.Filter("evtfilter_track_selection_cut==true");
    // Create the RDF for PSD charge analysis
    auto _fr_psd_charge_analysis = _fr_stk_analysis.Filter("evtfilter_psd_charge_measurement==true");
    // Create the RDF for STK charge analysis
    auto _fr_stk_charge_analysis = _fr_stk_analysis.Filter("evtfilter_stk_charge_measurement==true");

    std::cout << "\n\n**** Filter statistics ****\n";
    std::cout << "***************************\n";
    std::cout << "\nBGO analysis events: " << *(_fr_bgo_analysis.Count());
    std::cout << "\nSTK analysis events: " << *(_fr_stk_analysis.Count());
    std::cout << "\nPSD charge analysis events: " << *(_fr_psd_charge_analysis.Count());
    std::cout << "\nSTK charge analysis events: " << *(_fr_stk_charge_analysis.Count());
    std::cout << "\nPreselected events: " << *(_fr_preselected.Count());
    std::cout << "\n\n***************************";

    if (verbose) std::cout << "\n\nAnlysis running...\n\n";

    // Get binning
    std::vector<float> sumRms_binning;
    std::vector<float> flast_binning;
    auto flast13_binning = createLogBinning(1e-5, 2e-1, 1e+3);
    regularize_vars ? sumRms_binning = createLinearBinning(-20, 20, 1e+3) : sumRms_binning = createLogBinning(10, 2e+3, 1e+2);
    regularize_vars ? flast_binning = createLinearBinning(-10, 10, 1e+3) : flast_binning = createLogBinning(1e-5, 2e-1, 1e+3);
    auto sumRms_cosine_bins = createLogBinning(10, 3e+3, 1e+2);
    auto xtrl_bins = createLinearBinning(0, 150, 1e+2);
    auto cosine_bins = createLinearBinning(0, 1, 1e+2);
    auto energy_ratio_bins = createLinearBinning(-1, 0.2, 1e+3);
    auto flast_zoom_binning = createLogBinning(1e-5, 1e-2, 1e+3);
    auto flast_cosine_binning = createLogBinning(1e-5, 3e-1, 1e+3);
    
    // Extract BGO histos
    auto h_BGOrec_raw_energy = _fr_bgo_analysis.Define("raw_energy_gev", "energy * 0.001")
                                   .Histo1D<double, double>({"h_BGOrec_raw_energy", "BGO Raw energy", energy_nbins, &energy_binning[0]}, "raw_energy_gev", "simu_energy_w_corr");
    auto h_BGOrec_corr_energy = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Histo1D<double, double>({"h_BGOrec_corr_energy", "BGO Corr energy", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "simu_energy_w_corr");
    auto h_BGOrec_cosine = _fr_bgo_analysis.Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                               .Histo1D<double, double>({"h_BGOrec_cosine", "BGO Reco cosine", 100, 0, 1}, "bgorec_cosine", "simu_energy_w_corr");

    auto h_BGOrec_sumRms_bin_cosine2D_20_100 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                   .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                                   .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                   .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_20_100", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_sumRms_bin_cosine2D_100_250 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                    .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                                    .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                    .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_100_250", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_sumRms_bin_cosine2D_250_500 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                    .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                                    .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                    .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_250_500", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_sumRms_bin_cosine2D_500_1000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                     .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                                     .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                     .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_500_1000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_sumRms_bin_cosine2D_1000_3000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                      .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                                      .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                      .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_1000_3000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_sumRms_bin_cosine2D_3000_10000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                       .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                                       .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                       .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_3000_10000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_sumRms_bin_cosine2D_10000_20000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                                        .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                        .Histo2D<double, double, double>({"h_BGOrec_sumRms_bin_cosine2D_10000_20000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");

    auto h_BGOrec_sumRms_flast_20_100 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                            .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                            .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_20_100", "sumRms vs F_{last} correlation - 20 GeV - 100 GeV; sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_100_250 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                             .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                             .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_100_250", "sumRms vs F_{last} correlation - 100 GeV - 250 GeV;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_250_500 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                             .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                             .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_250_500", "sumRms vs F_{last} correlation - 250 GeV - 500 GeV;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_500_1000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                              .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                              .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_500_1000", "sumRms vs F_{last} correlation - 500 GeV - 1 TeV;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_1000_3000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                               .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                               .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_1000_3000", "sumRms vs F_{last} correlation - 1 TeV - 3 TeV;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_3000_10000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                                .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_3000_10000", "sumRms vs F_{last} correlation - 3 TeV - 10 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_10000_20000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                 .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                                 .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_10000_20000", "sumRms vs F_{last} correlation - 10 TeV - 20 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");

    auto h_BGOrec_sumRms_flast_13_20_100 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                               .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                               .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_20_100", "sumRms vs F_{13} correlation - 20 GeV - 100 GeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_13_100_250 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                                .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_100_250", "sumRms vs F_{13} correlation - 100 GeV - 250 GeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_13_250_500 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                                .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_250_500", "sumRms vs F_{13} correlation - 250 GeV - 500 GeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_13_500_1000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                 .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                                 .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_500_1000", "sumRms vs F_{13} correlation - 500 GeV - 1 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_13_1000_3000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                  .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                                  .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_1000_3000", "sumRms vs F_{13} correlation - 1 TeV - 3 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_13_3000_10000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                   .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                                   .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_3000_10000", "sumRms vs F_{13} correlation - 3 TeV - 10 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_sumRms_flast_13_10000_20000 = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                                                    .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                                    .Histo2D<double, double, double>({"h_BGOrec_sumRms_flast_13_10000_20000", "sumRms vs F_{13} correlation - 10 TeV - 20 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_last_layer = _fr_bgo_analysis.Histo1D<int, double>({"h_BGOrec_last_layer", "Last BGO layer; Last Energy Layer; entries", 14, 0, 14}, "lastBGOLayer", "simu_energy_w_corr");

    std::shared_ptr<TH1D> h_BGOrec_bar_energy = std::make_shared<TH1D>("h_BGOrec_bar_energy", "BGO Bar Energy; Bar Energy [MeV]", 1e+3, 0, 1e+4);
    std::shared_ptr<TH1D> h_BGOrec_bar_energy_zoom = std::make_shared<TH1D>("h_BGOrec_bar_energy_zoom", "BGO Bar Energy; Bar Energy [MeV]", 1e+3, 0, 1e+2);

    std::shared_ptr<TH1D> h_BGOrec_bar_energy_10MeV = std::make_shared<TH1D>("h_BGOrec_bar_energy_10MeV", "BGO Bar Energy - 10 MeV filter; Bar Energy [MeV]", 1e+3, 10, 1e+4);
    std::shared_ptr<TH1D> h_BGOrec_bar_energy_20MeV = std::make_shared<TH1D>("h_BGOrec_bar_energy_20MeV", "BGO Bar Energy - 20 MeV filter; Bar Energy [MeV]", 1e+3, 20, 1e+4);
    std::shared_ptr<TH1D> h_BGOrec_bar_energy_30MeV = std::make_shared<TH1D>("h_BGOrec_bar_energy_30MeV", "BGO Bar Energy - 30 MeV filter; Bar Energy [MeV]", 1e+3, 30, 1e+4);
    std::shared_ptr<TH1D> h_BGOrec_bar_energy_50MeV = std::make_shared<TH1D>("h_BGOrec_bar_energy_50MeV", "BGO Bar Energy - 50 MeV filter; Bar Energy [MeV]", 1e+3, 50, 1e+4);
    std::shared_ptr<TH1D> h_BGOrec_bar_energy_100MeV = std::make_shared<TH1D>("h_BGOrec_bar_energy_100MeV", "BGO Bar Energy - 100 MeV filter; Bar Energy [MeV]", 1e+3, 100, 1e+4);

    auto h_BGOrec_hits = _fr_bgo_analysis.Histo1D<int, double>({"h_BGOrec_hits", "BGO Hits", 100, 0, 400}, "nBGOentries", "simu_energy_w_corr");
    std::shared_ptr<TH1D> h_BGOrec_hits_10MeV = std::make_shared<TH1D>("h_BGOrec_hits_10MeV", "BGO Hits - 10 MeV filter", 1e+2, 0, 4e+2);
    std::shared_ptr<TH1D> h_BGOrec_hits_20MeV = std::make_shared<TH1D>("h_BGOrec_hits_20MeV", "BGO Hits - 20 MeV filter", 1e+2, 0, 4e+2);
    std::shared_ptr<TH1D> h_BGOrec_hits_30MeV = std::make_shared<TH1D>("h_BGOrec_hits_30MeV", "BGO Hits - 30 MeV filter", 1e+2, 0, 4e+2);
    std::shared_ptr<TH1D> h_BGOrec_hits_50MeV = std::make_shared<TH1D>("h_BGOrec_hits_50MeV", "BGO Hits - 50 MeV filter", 1e+2, 0, 4e+2);
    std::shared_ptr<TH1D> h_BGOrec_hits_100MeV = std::make_shared<TH1D>("h_BGOrec_hits_100MeV", "BGO Hits - 100 MeV filter", 1e+2, 0, 4e+2);

    _fr_bgo_analysis.Foreach([
        &h_BGOrec_bar_energy,
        &h_BGOrec_bar_energy_zoom, 
        &h_BGOrec_bar_energy_10MeV, 
        &h_BGOrec_bar_energy_20MeV, 
        &h_BGOrec_bar_energy_30MeV, 
        &h_BGOrec_bar_energy_50MeV, 
        &h_BGOrec_bar_energy_100MeV,
        &h_BGOrec_hits_10MeV,
        &h_BGOrec_hits_20MeV,
        &h_BGOrec_hits_30MeV,
        &h_BGOrec_hits_50MeV,
        &h_BGOrec_hits_100MeV]
        (std::vector<std::vector<double>> bar_energy, double energy_w)
        {
            unsigned int hits_10=0;
            unsigned int hits_20=0;
            unsigned int hits_30=0;
            unsigned int hits_50=0;
            unsigned int hits_100=0;

            for (auto& bgo_ly : bar_energy)
                for(auto& energy : bgo_ly)
                { 
                    h_BGOrec_bar_energy->Fill(energy, energy_w);
                    if (energy <= 100)
                        h_BGOrec_bar_energy_zoom->Fill(energy, energy_w);
                    if (energy>=10)
                    {
                        h_BGOrec_bar_energy_10MeV->Fill(energy, energy_w);
                        ++hits_10;
                    }
                    if (energy>=20)
                    {
                        h_BGOrec_bar_energy_20MeV->Fill(energy, energy_w);
                        ++hits_20;
                    }
                    if (energy>=30)
                    {
                        h_BGOrec_bar_energy_30MeV->Fill(energy, energy_w);
                        ++hits_30;
                    }
                    if (energy>=50)
                    {
                        h_BGOrec_bar_energy_50MeV->Fill(energy, energy_w);
                        ++hits_50;
                    }
                    if (energy>=100)
                    {
                        h_BGOrec_bar_energy_100MeV->Fill(energy, energy_w);
                        ++hits_100;
                    }
                }
            
            h_BGOrec_hits_10MeV->Fill(hits_10, energy_w);
            h_BGOrec_hits_20MeV->Fill(hits_20, energy_w);
            h_BGOrec_hits_30MeV->Fill(hits_30, energy_w);
            h_BGOrec_hits_50MeV->Fill(hits_50, energy_w);
            h_BGOrec_hits_100MeV->Fill(hits_100, energy_w);

        } , {"layerBarEnergy", "simu_energy_w_corr"});


    auto h_BGOrec_slopeX = _fr_bgo_analysis.Histo1D<double, double>({"h_BGOrec_slopeX", "BGOrec Slope X", 200, -10, 10}, "BGOrec_slopeX", "simu_energy_w_corr");
    auto h_BGOrec_slopeY = _fr_bgo_analysis.Histo1D<double, double>({"h_BGOrec_slopeY", "BGOrec Slope Y", 200, -10, 10}, "BGOrec_slopeY", "simu_energy_w_corr");
    auto h_BGOrec_interceptX = _fr_bgo_analysis.Histo1D<double, double>({"h_BGOrec_interceptX", "BGOrec Intercept X", 100, -500, 500}, "BGOrec_interceptX", "simu_energy_w_corr");
    auto h_BGOrec_interceptY = _fr_bgo_analysis.Histo1D<double, double>({"h_BGOrec_interceptY", "BGOrec Intercept Y", 100, -500, 500}, "BGOrec_interceptY", "simu_energy_w_corr");

    auto h_xtrl_energy_int = _fr_bgo_analysis.Histo1D<double>({"h_xtrl_energy_int", "XTRL", 200, 0, 150}, "xtrl", "simu_energy_w_corr");
    auto h_xtrl = _fr_bgo_analysis.Define("corr_energy_gev", "energy_corr * 0.001")
                      .Histo2D<double, double, double>({"h_xtrl", "XTRL", energy_nbins, &energy_binning[0], (int)xtrl_bins.size() - 1, &xtrl_bins[0]}, "corr_energy_gev", "xtrl", "simu_energy_w_corr");

    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_cosine_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_sumRms_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_sumRms_mean_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_sumRms_bin_weight(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_sumRms_sumRms_weight(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_sumRms_bin_cosine(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_sumRms_bin_cosine_cosine2D(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_sumRms_bin_cosine_cosine2D_log(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_sumRms_bin_cosine2D(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_sumRms_bin_cosine2D_log(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_sumRms_flast_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_sumRms_flast_13_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_max_energy_ratio(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ratio_last(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ratio_13(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ratio_last_cosine(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ratio_last_cosine_cosine2D(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ratio_last_cosine_cosine2D_log(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ratio_last_cosine2D(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ratio_last_cosine2D_log(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ratio_last_cosine2D_fdr(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ratio_last_cosine2D_fdr_log(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_last_layer_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_hits_bin(energy_nbins);

    std::vector<std::shared_ptr<TH1D>> h_BGOrec_hits_bin_10MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_hits_bin_20MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_hits_bin_30MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_hits_bin_50MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_hits_bin_100MeV(energy_nbins);

    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_xtrl_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_stkclusters_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_stkclusters_bgohits_bin(energy_nbins);
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH2D>>> h_stkclusters_bgohits_layer(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH2D>>(DAMPE_stk_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_rms_layer(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_energy_ratio_layer(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_energy_ratio_1R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_energy_ratio_2R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_energy_ratio_3R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_energy_ratio_5R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));

    double sumrms_bin_lvalue, sumrms_bin_rvalue;
    double flast_bin_lvalue, flast_bin_rvalue;
    double sumrms2D_bin_lvalue, sumrms2D_bin_rvalue;

    regularize_vars ? sumrms_bin_lvalue = -20 : sumrms_bin_lvalue = 0;
    regularize_vars ? sumrms_bin_rvalue = 20 : sumrms_bin_rvalue = 1e+3;
    regularize_vars ? sumrms2D_bin_lvalue = -20 : sumrms2D_bin_lvalue = 0;
    regularize_vars ? sumrms2D_bin_rvalue = 20 : sumrms2D_bin_rvalue = 2e+3;
    regularize_vars ? flast_bin_lvalue = -10 : flast_bin_lvalue = 0;
    regularize_vars ? flast_bin_rvalue = 10 : flast_bin_rvalue = 0.4;

    std::vector<std::shared_ptr<TH1D>> h_BGOrec_bar_energy_bin(energy_nbins);
    
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_bar_energy_bin_10_MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_bar_energy_bin_20_MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_bar_energy_bin_30_MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_bar_energy_bin_50_MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_bar_energy_bin_100_MeV(energy_nbins);

    std::vector<std::shared_ptr<TH2D>> h_BGOrec_shower_profile(energy_nbins);
    std::vector<std::shared_ptr<TH2D>> h_BGOrec_shower_profile_upto_09(energy_nbins);
    std::vector<std::shared_ptr<TH2D>> h_BGOrec_shower_profile_from_09(energy_nbins);

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        h_BGOrec_bar_energy_bin[bin_idx - 1] = std::make_shared<TH1D>((std::string("h_BGOrec_bar_energy_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 1e+3, 0, 1e+4);
        h_BGOrec_shower_profile[bin_idx - 1] = std::make_shared<TH2D>((std::string("h_BGOrec_shower_profile_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Shower Profile - bin ") + std::to_string(bin_idx)).c_str(), 13, 0, 13, 100, 0, 1);
        h_BGOrec_shower_profile_upto_09[bin_idx - 1] = std::make_shared<TH2D>((std::string("h_BGOrec_shower_profile_upto_09_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Shower Profile cos(#theta) < 0.9 - bin ") + std::to_string(bin_idx)).c_str(), 13, 0, 13, 100, 0, 1);
        h_BGOrec_shower_profile_from_09[bin_idx - 1] = std::make_shared<TH2D>((std::string("h_BGOrec_shower_profile_from_09_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Shower Profile cos(#theta) > 0.9 - bin ") + std::to_string(bin_idx)).c_str(), 13, 0, 13, 100, 0, 1);
    
        h_BGOrec_hits_bin_10MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_hits_10MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Hits - 10 MeV filter - bin ") + std::to_string(bin_idx)).c_str(), 1e+2, 0, 4e+2);
        h_BGOrec_hits_bin_20MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_hits_20MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Hits - 20 MeV filter - bin ") + std::to_string(bin_idx)).c_str(), 1e+2, 0, 4e+2);
        h_BGOrec_hits_bin_30MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_hits_30MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Hits - 30 MeV filter - bin ") + std::to_string(bin_idx)).c_str(), 1e+2, 0, 4e+2);
        h_BGOrec_hits_bin_50MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_hits_50MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Hits - 50 MeV filter - bin ") + std::to_string(bin_idx)).c_str(), 1e+2, 0, 4e+2);
        h_BGOrec_hits_bin_100MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_hits_100MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Hits - 100 MeV filter - bin ") + std::to_string(bin_idx)).c_str(), 1e+2, 0, 4e+2);

        h_BGOrec_bar_energy_bin_10_MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_bar_energy_10_MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - 10 MeV filter - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 1e+3, 10, 1e+4);
        h_BGOrec_bar_energy_bin_20_MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_bar_energy_20_MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - 20 MeV filter - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 1e+3, 20, 1e+4);
        h_BGOrec_bar_energy_bin_30_MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_bar_energy_30_MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - 30 MeV filter - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 1e+3, 30, 1e+4);
        h_BGOrec_bar_energy_bin_50_MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_bar_energy_50_MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - 50 MeV filter - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 1e+3, 50, 1e+4);
        h_BGOrec_bar_energy_bin_100_MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_bar_energy_100_MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - 100 MeV filter - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 1e+3, 100, 1e+4);
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
                                               .Histo1D<double, double>({(std::string("h_BGOrec_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Reco cosine - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "bgorec_cosine", "simu_energy_w_corr");
        h_BGOrec_sumRms_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                               .Histo1D<double, double>({(std::string("h_BGOrec_sumRms_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms - bin ") + std::to_string(bin_idx)).c_str(), 100, sumrms_bin_lvalue, sumrms_bin_rvalue}, sumRms_leaf.c_str(), "simu_energy_w_corr");
        h_BGOrec_sumRms_mean_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                    .Define("num_layers", [](int last_layer) -> int { return last_layer + 1; }, {"lastBGOLayer"})
                                                    .Define("sumRms_eff", "sumRms/(lastBGOLayer+1)")
                                                    .Histo1D<double, double>({(std::string("h_BGOrec_sumRms_mean_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms - bin ") + std::to_string(bin_idx)).c_str(), 80, 0, 100}, "sumRms_eff", "simu_energy_w_corr");
        h_BGOrec_sumRms_bin_weight[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                      .Define("sumRms_w", computeSumRmsWeight, {"eLayer", "rmsLayer", "energy"})
                                                      .Histo1D<double, double>({(std::string("h_BGOrec_sumRms_weight_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms weighted - bin ") + std::to_string(bin_idx)).c_str(), 50, 0, 100}, "sumRms_w", "simu_energy_w_corr");
        h_BGOrec_sumRms_sumRms_weight[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                         .Define("sumRms_w", computeSumRmsWeight, {"eLayer", "rmsLayer", "energy"})
                                                         .Define("sumRms_eff", "sumRms/(lastBGOLayer+1)")
                                                         .Histo2D<double, double, double>({(std::string("h_BGOrec_sumRms_sumRms_weight_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms_{mean} vs sumRms_{weighted} - bin ") + std::to_string(bin_idx) + std::string("; sumRms_{weighted} [mm]; sumRms_{mean} [mm]")).c_str(), 50, 0, 100, 50, 0, 100}, "sumRms_w", "sumRms_eff", "simu_energy_w_corr");
        h_BGOrec_sumRms_bin_cosine[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                      .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                      .Define("sumRms_cosine", "sumRms/bgorec_cosine")
                                                      .Histo1D<double, double>({(std::string("h_BGOrec_sumRms_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; sumRms/cos(#theta); counts")).c_str(), 100, 0, 3000}, "sumRms_cosine", "simu_energy_w_corr");
        h_BGOrec_sumRms_bin_cosine_cosine2D[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                               .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                               .Define("sumRms_cosine", "sumRms/bgorec_cosine")
                                                               .Histo2D<double, double, double>({(std::string("h_BGOrec_sumRms_cosine_cosine2D_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); sumRms/cos(#theta)")).c_str(), 100, 0, 1, 100, 0, 3000}, "bgorec_cosine", "sumRms_cosine", "simu_energy_w_corr");
        h_BGOrec_sumRms_bin_cosine_cosine2D_log[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                                   .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                   .Define("sumRms_cosine", "sumRms/bgorec_cosine")
                                                                   .Histo2D<double, double, double>({(std::string("h_BGOrec_sumRms_cosine_cosine2D_log_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); sumRms/cos(#theta)")).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_cosine_bins.size() - 1, &sumRms_cosine_bins[0]}, "bgorec_cosine", "sumRms_cosine", "simu_energy_w_corr");
        h_BGOrec_sumRms_bin_cosine2D[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                        .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                        .Histo2D<double, double>({(std::string("h_BGOrec_sumRms_cosine2D_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); sumRms [mm]")).c_str(), 100, 0, 1, 100, sumrms2D_bin_lvalue, sumrms2D_bin_rvalue}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
        if (!regularize_vars)
            h_BGOrec_sumRms_bin_cosine2D_log[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                .Histo2D<double, double>({(std::string("h_BGOrec_sumRms_cosine2D_log_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); sumRms [mm]")).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
        h_BGOrec_sumRms_flast_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                     .Histo2D<double, double, double>({(std::string("h_BGOrec_sumRms_flast_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms vs F_{last} - bin ") + std::to_string(bin_idx) + std::string("; sumRms [mm]; F_{last}")).c_str(), (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
        h_BGOrec_sumRms_flast_13_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                        .Histo2D<double, double, double>({(std::string("h_BGOrec_sumRms_flast_13_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms vs F_{13} - bin ") + std::to_string(bin_idx) + std::string("; sumRms [mm]; F_{13}")).c_str(), (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
        h_BGOrec_max_energy_ratio[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                     .Define("max_energy_ratio", GetMaxEnergyRatio, {"fracLayer"})
                                                     .Histo1D<double, double>({(std::string("h_BGOrec_max_energy_ratio_bin_") + std::to_string(bin_idx)).c_str(), (std::string("max energy ratio - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "max_energy_ratio", "simu_energy_w_corr");
        h_BGOrec_ratio_last[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                               .Histo1D<double, double>({(std::string("h_BGOrec_ratio_last_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer - bin ") + std::to_string(bin_idx)).c_str(), 100, flast_bin_lvalue, flast_bin_rvalue}, fracLast_leaf, "simu_energy_w_corr");
        h_BGOrec_ratio_13[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                             .Define("energy_ratio_13", "fracLayer[13]")
                                             .Histo1D<double, double>({(std::string("h_BGOrec_ratio_13_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 13th Layer - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 0.01}, "energy_ratio_13", "simu_energy_w_corr");
        h_BGOrec_ratio_last_cosine[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                      .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                      .Define("frac_last_cosine", "fracLast/bgorec_cosine")
                                                      .Histo1D<double, double>({(std::string("h_BGOrec_ratio_last_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("F_{last}/cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; F_{last}/cos(#theta); entries")).c_str(), 1000, 0, 0.3}, "frac_last_cosine", "simu_energy_w_corr");
        h_BGOrec_ratio_last_cosine_cosine2D[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                               .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                               .Define("frac_last_cosine", "fracLast/bgorec_cosine")
                                                               .Histo2D<double, double, double>({(std::string("h_BGOrec_ratio_last_cosine_cosine2D_bin_") + std::to_string(bin_idx)).c_str(), (std::string("F_{last}/cos(#theta) vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}/cos(#theta)")).c_str(), 100, 0, 1, 1000, 0, 0.3}, "bgorec_cosine", "frac_last_cosine", "simu_energy_w_corr");
        h_BGOrec_ratio_last_cosine_cosine2D_log[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                                   .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                   .Define("frac_last_cosine", "fracLast/bgorec_cosine")
                                                                   .Histo2D<double, double, double>({(std::string("h_BGOrec_ratio_last_cosine_cosine2D_log_bin_") + std::to_string(bin_idx)).c_str(), (std::string("F_{last}/cos(#theta) vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}/cos(#theta)")).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)flast_cosine_binning.size() - 1, &flast_cosine_binning[0]}, "bgorec_cosine", "frac_last_cosine", "simu_energy_w_corr");
        h_BGOrec_ratio_last_cosine2D[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                        .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                        .Histo2D<double, double, double>({(std::string("h_BGOrec_ratio_last_cosine2D_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}")).c_str(), 100, 0, 1, 100, flast_bin_lvalue, flast_bin_rvalue}, "bgorec_cosine", fracLast_leaf, "simu_energy_w_corr");
        if (!regularize_vars)
        {
            h_BGOrec_ratio_last_cosine2D_log[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                .Histo2D<double, double, double>({(std::string("h_BGOrec_ratio_last_cosine2D_log_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}")).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)flast_zoom_binning.size() - 1, &flast_zoom_binning[0]}, "bgorec_cosine", fracLast_leaf, "simu_energy_w_corr");
            h_BGOrec_ratio_last_cosine2D_fdr[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                .Histo2D<double, double, double>({(std::string("h_BGOrec_ratio_last_cosine2D_fdr_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}")).c_str(), 100, 0, 1, 1000, 0, 0.2}, "bgorec_cosine", fracLast_leaf, "simu_energy_w_corr");
            h_BGOrec_ratio_last_cosine2D_fdr_log[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                                    .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                    .Histo2D<double, double, double>({(std::string("h_BGOrec_ratio_last_cosine2D_fdr_log_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}")).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "bgorec_cosine", fracLast_leaf, "simu_energy_w_corr");
        }
        h_BGOrec_last_layer_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                   .Histo1D<int, double>({(std::string("h_BGOrec_last_layer_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Last BGO layer - bin ") + std::to_string(bin_idx) + std::string("; Last Energy Layer; entries")).c_str(), 14, 0, 14}, "lastBGOLayer", "simu_energy_w_corr");
        h_BGOrec_hits_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                         .Histo1D<int, double>({(std::string("h_BGOrec_hits_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO hits - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 400}, "nBGOentries", "simu_energy_w_corr");
        h_xtrl_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                      .Histo1D<double, double>({(std::string("h_xtrl_bin_") + std::to_string(bin_idx)).c_str(), (std::string("XTRL - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 150}, "xtrl", "simu_energy_w_corr");
        h_stkclusters_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                        .Define("stk_clusters",
                                                                 [](const std::vector<int> plane_clusters) -> int { return std::accumulate(plane_clusters.begin(), plane_clusters.end(), 0); }, {"STK_plane_clusters"})
                                             .Histo1D<int, double>({(std::string("h_stkclusters_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Total STK clusters - bin ") + std::to_string(bin_idx) + std::string(" ; total STK clusters; entries")).c_str(), 100, 0, 1500}, "stk_clusters", "simu_energy_w_corr");
        h_stkclusters_bgohits_bin[bin_idx - 1] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                                .Define("stk_clusters",
                                                                         [](const std::vector<int> plane_clusters) -> int { return std::accumulate(plane_clusters.begin(), plane_clusters.end(), 0); }, {"STK_plane_clusters"})
                                                     .Histo2D<int, int, double>({(std::string("h_stkclusters_bgohist_bin_") + std::to_string(bin_idx)).c_str(), (std::string("STK clusters vs BGO hits - bin ") + std::to_string(bin_idx) + std::string(" ; STK clusters; BGO hits")).c_str(), 100, 0, 500, 80, 0, 350}, "stk_clusters", "nBGOentries", "simu_energy_w_corr");

        _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
            .Foreach([
                &h_BGOrec_bar_energy_bin,
                &h_BGOrec_bar_energy_bin_10_MeV,
                &h_BGOrec_bar_energy_bin_20_MeV,
                &h_BGOrec_bar_energy_bin_30_MeV,
                &h_BGOrec_bar_energy_bin_50_MeV,
                &h_BGOrec_bar_energy_bin_100_MeV,
                &h_BGOrec_hits_bin_10MeV,
                &h_BGOrec_hits_bin_20MeV,
                &h_BGOrec_hits_bin_30MeV,
                &h_BGOrec_hits_bin_50MeV,
                &h_BGOrec_hits_bin_100MeV,
                bin_idx]
                (std::vector<std::vector<double>> bar_energy, double energy_w) 
                { 
                    unsigned int hits_10=0;
                    unsigned int hits_20=0;
                    unsigned int hits_30=0;
                    unsigned int hits_50=0;
                    unsigned int hits_100=0;

                    for (auto& bgo_ly : bar_energy)
                        for(auto& energy : bgo_ly)
                        {
                            h_BGOrec_bar_energy_bin[bin_idx -1]->Fill(energy, energy_w); 

                            if (energy>=10)
                            {
                                h_BGOrec_bar_energy_bin_10_MeV[bin_idx -1]->Fill(energy, energy_w);
                                ++hits_10;
                            }
                            if (energy>=20)
                            {
                                h_BGOrec_bar_energy_bin_20_MeV[bin_idx -1]->Fill(energy, energy_w);
                                ++hits_20;
                            }
                            if (energy>=30)
                            {
                                h_BGOrec_bar_energy_bin_30_MeV[bin_idx -1]->Fill(energy, energy_w);
                                ++hits_30;
                            }
                            if (energy>=50)
                            {
                                h_BGOrec_bar_energy_bin_50_MeV[bin_idx -1]->Fill(energy, energy_w);
                                ++hits_50;
                            }
                            if (energy>=100)
                            {
                                h_BGOrec_bar_energy_bin_100_MeV[bin_idx -1]->Fill(energy, energy_w);
                                ++hits_50;
                            }
                        }
                    
                    h_BGOrec_hits_bin_10MeV[bin_idx -1]->Fill(hits_10, energy_w);
                    h_BGOrec_hits_bin_20MeV[bin_idx -1]->Fill(hits_20, energy_w);
                    h_BGOrec_hits_bin_30MeV[bin_idx -1]->Fill(hits_30, energy_w);
                    h_BGOrec_hits_bin_50MeV[bin_idx -1]->Fill(hits_50, energy_w);
                    h_BGOrec_hits_bin_100MeV[bin_idx -1]->Fill(hits_100, energy_w);

                }, {"layerBarEnergy", "simu_energy_w_corr"});

        _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
            .Foreach([&h_BGOrec_shower_profile, bin_idx](std::vector<double> &energy_frac, double energy_w) { 
                for (int lidx=0; lidx<DAMPE_bgo_nLayers; ++lidx)
                   h_BGOrec_shower_profile[bin_idx - 1]->Fill(lidx, energy_frac[lidx], energy_w); }, {"fracLayer", "simu_energy_w_corr"});

        _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
            .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
            .Filter("bgorec_cosine<0.9")
            .Foreach([&h_BGOrec_shower_profile_upto_09, bin_idx](std::vector<double> &energy_frac, double energy_w) { 
                for (int lidx=0; lidx<DAMPE_bgo_nLayers; ++lidx)
                    h_BGOrec_shower_profile_upto_09[bin_idx - 1]->Fill(lidx, energy_frac[lidx], energy_w); }, {"fracLayer", "simu_energy_w_corr"});

        _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
            .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
            .Filter("bgorec_cosine>0.9")
            .Foreach([&h_BGOrec_shower_profile_from_09, bin_idx](std::vector<double> &energy_frac, double energy_w) { 
                for (int lidx=0; lidx<DAMPE_bgo_nLayers; ++lidx)
                    h_BGOrec_shower_profile_from_09[bin_idx - 1]->Fill(lidx, energy_frac[lidx], energy_w); }, {"fracLayer", "simu_energy_w_corr"});

        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto GetLayerComponent = [=](std::vector<double> &value_layer) -> double { return value_layer[ly]; };
            auto GetStkLayerComponent = [=](std::vector<int> &stk_clusters) -> int { return stk_clusters[ly]; };

            if (ly < DAMPE_stk_nLayers)
                h_stkclusters_bgohits_layer[bin_idx -1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                                .Define("stk_layer_clusters", GetStkLayerComponent, {"STK_plane_clusters"})
                                                                .Histo2D<int, int, double>({(std::string("h_stkclusters_bgohits_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("STK clusters vs BGO hits - bin ") + std::to_string(bin_idx) + std::string(" ; STK clusters; BGO hits")).c_str(), 100, 0, 500, 80, 0, 350}, "stk_layer_clusters", "nBGOentries", "simu_energy_w_corr");
            h_BGOrec_rms_layer[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                      .Define("rms_layer", GetLayerComponent, {"rmsLayer"})
                                                      .Histo1D<double, double>({(std::string("h_BGOrec_rms_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("RMS layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 500}, "rms_layer", "simu_energy_w_corr");
            h_BGOrec_energy_ratio_layer[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                               .Define("energy_ratio_layer", GetLayerComponent, {"fracLayer"})
                                                               .Histo1D<double, double>({(std::string("h_BGOrec_energy_ratio_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_layer", "simu_energy_w_corr");
            h_BGOrec_energy_ratio_1R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                            .Define("energy_1R_layer", GetLayerComponent, {"energy_1R_radius"})
                                                            .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                            .Define("energy_ratio_1R_layer", "energy_1R_layer/energy_layer")
                                                            .Histo1D<double, double>({(std::string("h_BGOrec_energy_ratio_1R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 1MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_1R_layer", "simu_energy_w_corr");
            h_BGOrec_energy_ratio_2R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                            .Define("energy_2R_layer", GetLayerComponent, {"energy_2R_radius"})
                                                            .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                            .Define("energy_ratio_2R_layer", "energy_2R_layer/energy_layer")
                                                            .Histo1D<double, double>({(std::string("h_BGOrec_energy_ratio_2R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 2MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_2R_layer", "simu_energy_w_corr");
            h_BGOrec_energy_ratio_3R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                            .Define("energy_3R_layer", GetLayerComponent, {"energy_3R_radius"})
                                                            .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                            .Define("energy_ratio_3R_layer", "energy_3R_layer/energy_layer")
                                                            .Histo1D<double, double>({(std::string("h_BGOrec_energy_ratio_3R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 3MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_3R_layer", "simu_energy_w_corr");
            h_BGOrec_energy_ratio_5R[bin_idx - 1][ly] = _fr_bgo_analysis.Filter(bin_filter, {"energy_bin"})
                                                            .Define("energy_5R_layer", GetLayerComponent, {"energy_5R_radius"})
                                                            .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                            .Define("energy_ratio_5R_layer", "energy_5R_layer/energy_layer")
                                                            .Histo1D<double, double>({(std::string("h_BGOrec_energy_ratio_5R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 5MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_5R_layer", "simu_energy_w_corr");
        }
    }

    ROOT::RDF::RResultPtr<TH1D> h_simu_energy;
    ROOT::RDF::RResultPtr<TH1D> h_simu_energy_w;
    ROOT::RDF::RResultPtr<TH1D> h_energy_diff;
    ROOT::RDF::RResultPtr<TH1D> h_energy_diff_corr;
    ROOT::RDF::RResultPtr<TH2D> h_energy_diff2D;
    ROOT::RDF::RResultPtr<TH2D> h_energy_diff2D_corr;
    ROOT::RDF::RResultPtr<TH2D> h_energy_unfold;
    ROOT::RDF::RResultPtr<TH2D> h_energy_unfold_corr;
    ROOT::RDF::RResultPtr<TH1D> h_simu_slopeX;
    ROOT::RDF::RResultPtr<TH1D> h_simu_slopeY;
    ROOT::RDF::RResultPtr<TH1D> h_simu_interceptX;
    ROOT::RDF::RResultPtr<TH1D> h_simu_interceptY;
    ROOT::RDF::RResultPtr<TH1D> h_simu_position_x;
    ROOT::RDF::RResultPtr<TH1D> h_simu_position_y;
    ROOT::RDF::RResultPtr<TH1D> h_simu_position_z;
    ROOT::RDF::RResultPtr<TH3D> h_simu_position;
    ROOT::RDF::RResultPtr<TH1D> h_simu_flux_w;
    ROOT::RDF::RResultPtr<TH1D> h_simu_w;
    ROOT::RDF::RResultPtr<TH1D> h_simu_particle;
    ROOT::RDF::RResultPtr<TH1D> h_simu_theta;
    ROOT::RDF::RResultPtr<TH1D> h_simu_phi;
    ROOT::RDF::RResultPtr<TH1D> h_simu_charge;
    ROOT::RDF::RResultPtr<TH1D> h_simu_cosx;
    ROOT::RDF::RResultPtr<TH1D> h_simu_cosy;
    ROOT::RDF::RResultPtr<TH1D> h_simu_cosz;
    ROOT::RDF::RResultPtr<TH1D> h_simu_zenith;
    ROOT::RDF::RResultPtr<TH1D> h_simu_azimuth;

    std::unique_ptr<TH1D> h_simu_thruthtrajectory_start_x;
    std::unique_ptr<TH1D> h_simu_thruthtrajectory_start_y;
    std::unique_ptr<TH1D> h_simu_thruthtrajectory_start_z;
    std::unique_ptr<TH1D> h_simu_thruthtrajectory_stop_x;
    std::unique_ptr<TH1D> h_simu_thruthtrajectory_stop_y;
    std::unique_ptr<TH1D> h_simu_thruthtrajectory_stop_z;
    std::unique_ptr<TH1D> h_simu_thruthtrajectory_trackID;
    std::unique_ptr<TH1D> h_simu_thruthtrajectory_parentID;
    std::unique_ptr<TH1D> h_simu_thruthtrajectory_pdgID;
    std::unique_ptr<TH1D> h_simu_thruthtrajectory_stop_index;

    if (_mc)
    {
        // Extract Simu histos
        h_simu_energy = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                                .Histo1D<double, double>({"h_simu_energy", "Simu energy", energy_nbins, &energy_binning[0]}, "simu_energy_gev", "simu_energy_w_corr");
        h_simu_energy_w = _fr_bgo_analysis.Histo1D<double>({"h_simu_energy_w", "Simu energy weight", 100, 0, 1}, "simu_energy_w_corr");
        h_energy_diff = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                                .Define("raw_energy_gev", "energy * 0.001")
                                .Define("energy_diff", "simu_energy_gev - raw_energy_gev")
                                .Histo1D<double, double>({"h_energy_diff", "Simu vs Raw Reco BGO energy: Simu Energy - Raw Energy (GeV); counts", 100, 0, 100}, "energy_diff", "simu_energy_w_corr");
        h_energy_diff_corr = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("energy_diff", "simu_energy_gev - corr_energy_gev")
                                    .Histo1D<double, double>({"h_energy_diff_corr", "Simu vs Corrected Reco BGO energy: Simu Energy - Corrected Energy (GeV); counts", 100, -100, 100}, "energy_diff", "simu_energy_w_corr");
        h_energy_diff2D = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                                .Define("energy_ratio", "(energy - simu_energy)/simu_energy")
                                .Histo2D<double, double, double>({"h_energy_diff2D", "Energy Ratio; Real Energy (GeV); (Raw - Simu)/Simu", energy_nbins, &energy_binning[0], (int)energy_ratio_bins.size() - 1, &(energy_ratio_bins[0])}, "simu_energy_gev", "energy_ratio", "simu_energy_w_corr");
        h_energy_diff2D_corr = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                                        .Define("energy_ratio", "(energy_corr-simu_energy)/simu_energy")
                                        .Histo2D<double, double, double>({"h_energy_diff2D_corr", "Energy Ratio; Real Energy (GeV); (Corr - Simu)/Simu", energy_nbins, &energy_binning[0], (int)energy_ratio_bins.size() - 1, &(energy_ratio_bins[0])}, "simu_energy_gev", "energy_ratio", "simu_energy_w_corr");
        h_energy_unfold = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                                .Define("raw_energy_gev", "energy * 0.001")
                                .Histo2D<double, double, double>({"h_energy_unfold", "Energy Unfolding Matrix; Real Energy (GeV); Raw Energy (GeV)", energy_nbins, &energy_binning[0], energy_nbins, &energy_binning[0]}, "simu_energy_gev", "raw_energy_gev", "simu_energy_w_corr");
        h_energy_unfold_corr = _fr_bgo_analysis.Define("simu_energy_gev", "simu_energy * 0.001")
                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Histo2D<double, double, double>({"h_energy_unfold_corr", "Energy Unfolding Matrix; Real Energy (GeV); Corr Energy (GeV)", energy_nbins, &energy_binning[0], energy_nbins, &energy_binning[0]}, "simu_energy_gev", "corr_energy_gev", "simu_energy_w_corr");
        h_simu_slopeX = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_slopeX", "Simu Slope X", 100, -10, 10}, "simu_slope_x", "simu_energy_w_corr");
        h_simu_slopeY = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_slopeY", "Simu Slope Y", 100, -10, 10}, "simu_slope_y", "simu_energy_w_corr");
        h_simu_interceptX = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_interceptX", "Simu Intercept X", 100, -500, 500}, "simu_intercept_x", "simu_energy_w_corr");
        h_simu_interceptY = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_interceptY", "Simu Intercept Y", 100, -500, 500}, "simu_intercept_y", "simu_energy_w_corr");
        h_simu_position_x = _fr_bgo_analysis.Define("simu_position_comp", "simu_position.fX")
                                    .Histo1D<double, double>({"h_simu_position_x", "Simu Position X; X [mm]", 100, -1500, 1500}, "simu_position_comp", "simu_energy_w_corr");
        h_simu_position_y = _fr_bgo_analysis.Define("simu_position_comp", "simu_position.fY")
                                    .Histo1D<double, double>({"h_simu_position_y", "Simu Position Y; Y [mm]", 100, -1500, 1500}, "simu_position_comp", "simu_energy_w_corr");
        h_simu_position_z = _fr_bgo_analysis.Define("simu_position_comp", "simu_position.fZ")
                                    .Histo1D<double, double>({"h_simu_position_z", "Simu Position Z; Z [mm]", 100, -1500, 1500}, "simu_position_comp", "simu_energy_w_corr");
        h_simu_position = _fr_bgo_analysis.Define("simu_position_comp_x", "simu_position.fX")
                                .Define("simu_position_comp_y", "simu_position.fY")
                                .Define("simu_position_comp_z", "simu_position.fZ")
                                .Histo3D<double, double, double, double>({"h_simu_position", "Simu Position; Y [mm]; X [mm]; Z [mm]", 100, -1500, 1500, 100, -1500, 1500, 100, -1500, 1500}, "simu_position_comp_y", "simu_position_comp_x", "simu_position_comp_z", "simu_energy_w_corr");
        h_simu_flux_w = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_flux_w", "Simu F_{w}; F_{w}", 100, 0, 2}, "simu_flux_w", "simu_energy_w_corr");
        h_simu_w = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_w", "Simu w; w", 100, 0, 2}, "simu_w", "simu_energy_w_corr");
        h_simu_particle = _fr_bgo_analysis.Histo1D<int, double>({"h_simu_particle", "Simu particle", 3, 0, 3}, "simu_n_particle", "simu_energy_w_corr");
        h_simu_theta = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_theta", "Simu #theta; #theta [deg]", 100, 0, 180}, "simu_theta", "simu_energy_w_corr");
        h_simu_phi = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_phi", "Simu #phi; #phi [deg]", 200, -180, 180}, "simu_phi", "simu_energy_w_corr");
        h_simu_charge = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_charge", "Simu charge;", 100, -2, 2}, "simu_charge", "simu_energy_w_corr");
        h_simu_cosx = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_cosx", "Simu #cos(x);", 100, -1, 1}, "simu_cos_x", "simu_energy_w_corr");
        h_simu_cosy = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_cosy", "Simu #cos(y);", 100, -1, 1}, "simu_cos_y", "simu_energy_w_corr");
        h_simu_cosz = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_cosz", "Simu #cos(z);", 100, -1, 1}, "simu_cos_z", "simu_energy_w_corr");
        h_simu_zenith = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_zenith", "Simu zenith; Zenith [deg]", 100, 0, 90}, "simu_zenith", "simu_energy_w_corr");
        h_simu_azimuth = _fr_bgo_analysis.Histo1D<double, double>({"h_simu_azimuth", "Simu azimuth; Azimuth [deg]", 200, -180, 180}, "simu_azimuth", "simu_energy_w_corr");

        h_simu_thruthtrajectory_start_x = std::make_unique<TH1D>("h_simu_thruthtrajectory_start_x", "Truth Trajectory Start X; X [mm]", 100, -1500, 1500);
        h_simu_thruthtrajectory_start_y = std::make_unique<TH1D>("h_simu_thruthtrajectory_start_y", "Truth Trajectory Start Y; Y [mm]", 100, -1500, 1500);
        h_simu_thruthtrajectory_start_z = std::make_unique<TH1D>("h_simu_thruthtrajectory_start_z", "Truth Trajectory Start Z; Z [mm]", 100, -1500, 1500);
        h_simu_thruthtrajectory_stop_x = std::make_unique<TH1D>("h_simu_thruthtrajectory_stop_x", "Truth Trajectory Stop X; X [mm]", 100, -3000, 3000);
        h_simu_thruthtrajectory_stop_y = std::make_unique<TH1D>("h_simu_thruthtrajectory_stop_y", "Truth Trajectory Stop Y; Y [mm]", 100, -3000, 3000);
        h_simu_thruthtrajectory_stop_z = std::make_unique<TH1D>("h_simu_thruthtrajectory_stop_z", "Truth Trajectory Stop Z; Z [mm]", 100, -3000, 3000);
        h_simu_thruthtrajectory_trackID = std::make_unique<TH1D>("h_simu_thruthtrajectory_trackID", "Truth Trajectory Track ID", 100, 0, 200);
        h_simu_thruthtrajectory_parentID = std::make_unique<TH1D>("h_simu_thruthtrajectory_parentID", "Truth Trajectory Parent ID", 2, 0, 1);
        h_simu_thruthtrajectory_pdgID = std::make_unique<TH1D>("h_simu_thruthtrajectory_pdgID", "Truth Trajectory PDG ID", 200, -2000, 3000);
        h_simu_thruthtrajectory_stop_index = std::make_unique<TH1D>("h_simu_thruthtrajectory_stop_index", "Truth Trajectory Stop Index", 100, -10, 10);

        _fr_bgo_analysis.Foreach([&h_simu_thruthtrajectory_start_x](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_thruthtrajectory_start_x->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_start_x", "simu_energy_w_corr"});
        _fr_bgo_analysis.Foreach([&h_simu_thruthtrajectory_start_y](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_thruthtrajectory_start_y->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_start_y", "simu_energy_w_corr"});
        _fr_bgo_analysis.Foreach([&h_simu_thruthtrajectory_start_z](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_thruthtrajectory_start_z->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_start_z", "simu_energy_w_corr"});
        _fr_bgo_analysis.Foreach([&h_simu_thruthtrajectory_stop_x](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_thruthtrajectory_stop_x->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_stop_x", "simu_energy_w_corr"});
        _fr_bgo_analysis.Foreach([&h_simu_thruthtrajectory_stop_y](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_thruthtrajectory_stop_y->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_stop_y", "simu_energy_w_corr"});
        _fr_bgo_analysis.Foreach([&h_simu_thruthtrajectory_stop_z](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_thruthtrajectory_stop_z->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_stop_z", "simu_energy_w_corr"});
        _fr_bgo_analysis.Foreach([&h_simu_thruthtrajectory_trackID](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_thruthtrajectory_trackID->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_trackID", "simu_energy_w_corr"});
        _fr_bgo_analysis.Foreach([&h_simu_thruthtrajectory_parentID](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_thruthtrajectory_parentID->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_parentID", "simu_energy_w_corr"});
        _fr_bgo_analysis.Foreach([&h_simu_thruthtrajectory_pdgID](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_thruthtrajectory_pdgID->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_PDG", "simu_energy_w_corr"});
        _fr_bgo_analysis.Foreach([&h_simu_thruthtrajectory_stop_index](std::vector<double> thruthtrajectory, double energy_w) { 
            for (auto& _elm : thruthtrajectory)
                h_simu_thruthtrajectory_stop_index->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_stop_index", "simu_energy_w_corr"});
    }

    // Extract STK histos
    auto h_stk_cosine = _fr_stk_analysis.Histo1D({"h_stk_cosine", "h_stk_cosine", 100, 0, 1}, "STK_bestTrack_costheta");
    auto h_stk_slopeX = _fr_stk_analysis.Histo1D<double, double>({"h_stk_slopeX", "Simu Slope X", 200, -10, 10}, "STK_bestTrack_slopeX", "simu_energy_w_corr");
    auto h_stk_slopeY = _fr_stk_analysis.Histo1D<double, double>({"h_stk_slopeY", "Simu Slope Y", 200, -10, 10}, "STK_bestTrack_slopeY", "simu_energy_w_corr");
    auto h_stk_interceptX = _fr_stk_analysis.Histo1D<double, double>({"h_stk_interceptX", "Simu Intercept X", 100, -500, 500}, "STK_bestTrack_interceptX", "simu_energy_w_corr");
    auto h_stk_interceptY = _fr_stk_analysis.Histo1D<double, double>({"h_stk_interceptY", "Simu Intercept Y", 100, -500, 500}, "STK_bestTrack_interceptY", "simu_energy_w_corr");

    // Extract PSD charge histos
    auto h_psd_chargeX = _fr_psd_charge_analysis.Histo1D<double, double>({"h_psd_chargeX", "PSD Charge X", 100, 0, 20}, "PSD_chargeX", "simu_energy_w_corr");
    auto h_psd_chargeY = _fr_psd_charge_analysis.Histo1D<double, double>({"h_psd_chargeY", "PSD Charge Y", 100, 0, 20}, "PSD_chargeY", "simu_energy_w_corr");
    auto h_psd_charge = _fr_psd_charge_analysis.Define("psd_charge", [](double psd_chargeX, double psd_chargeY) { return 0.5 * (psd_chargeX + psd_chargeY); }, {"PSD_chargeX", "PSD_chargeY"})
                            .Histo1D<double, double>({"h_psd_charge", "PSD Charge", 100, 0, 20}, "psd_charge", "simu_energy_w_corr");
    auto h_psd_charge2D = _fr_psd_charge_analysis.Histo2D<double, double, double>({"h_psd_charge2D", "PSD Charge", 100, 0, 20, 100, 0, 20}, "PSD_chargeX", "PSD_chargeY", "simu_energy_w_corr");

    auto h_psd_selected_chargeX = _fr_psd_charge_analysis.Filter("evtfilter_psd_charge_cut==true")
                                      .Histo1D<double, double>({"h_psd_selected_chargeX", "PSD Charge X", 100, 0, 20}, "PSD_chargeX", "simu_energy_w_corr");
    auto h_psd_selected_chargeY = _fr_psd_charge_analysis.Filter("evtfilter_psd_charge_cut==true")
                                      .Histo1D<double, double>({"h_psd_selected_chargeY", "PSD Charge Y", 100, 0, 20}, "PSD_chargeY", "simu_energy_w_corr");
    auto h_psd_selected_charge = _fr_psd_charge_analysis.Filter("evtfilter_psd_charge_cut==true")
                                     .Define("psd_charge", [](double psd_chargeX, double psd_chargeY) { return 0.5 * (psd_chargeX + psd_chargeY); }, {"PSD_chargeX", "PSD_chargeY"})
                                     .Histo1D<double, double>({"h_psd_selected_charge", "PSD Charge", 100, 0, 20}, "psd_charge", "simu_energy_w_corr");
    auto h_psd_selected_charge2D = _fr_psd_charge_analysis.Filter("evtfilter_psd_charge_cut==true")
                                       .Histo2D<double, double, double>({"h_psd_selected_charge2D", "PSD Charge", 100, 0, 20, 100, 0, 20}, "PSD_chargeX", "PSD_chargeY", "simu_energy_w_corr");

    // Extract STK charge histos
    auto h_stk_chargeX = _fr_stk_charge_analysis.Histo1D<double, double>({"h_stk_chargeX", "STK Charge X", 500, 0, 100}, "STK_chargeX", "simu_energy_w_corr");
    auto h_stk_chargeY = _fr_stk_charge_analysis.Histo1D<double, double>({"h_stk_chargeY", "STK Charge Y", 500, 0, 100}, "STK_chargeY", "simu_energy_w_corr");
    auto h_stk_charge = _fr_stk_charge_analysis.Define("stk_charge", [](double stk_chargeX, double stk_chargeY) { return 0.5 * (stk_chargeX + stk_chargeY); }, {"STK_chargeX", "STK_chargeY"})
                            .Histo1D<double, double>({"h_stk_charge", "STK Charge", 100, 0, 20}, "stk_charge", "simu_energy_w_corr");
    auto h_stk_charge2D = _fr_stk_charge_analysis.Histo2D<double, double, double>({"h_stk_charge2D", "STK Charge", 500, 0, 100, 500, 0, 100}, "STK_chargeX", "STK_chargeY", "simu_energy_w_corr");

    auto h_stk_selected_chargeX = _fr_stk_charge_analysis.Filter("evtfilter_stk_charge_cut==true")
                                      .Histo1D<double, double>({"h_stk_selected_chargeX", "STK Charge X", 500, 0, 100}, "STK_chargeX", "simu_energy_w_corr");
    auto h_stk_selected_chargeY = _fr_stk_charge_analysis.Filter("evtfilter_stk_charge_cut==true")
                                      .Histo1D<double, double>({"h_stk_selected_chargeY", "STK Charge Y", 500, 0, 100}, "STK_chargeY", "simu_energy_w_corr");
    auto h_stk_selected_charge = _fr_stk_charge_analysis.Filter("evtfilter_stk_charge_cut==true")
                                     .Define("stk_charge", [](double stk_chargeX, double stk_chargeY) { return 0.5 * (stk_chargeX + stk_chargeY); }, {"STK_chargeX", "STK_chargeY"})
                                     .Histo1D<double, double>({"h_stk_selected_charge", "STK Charge", 500, 0, 100}, "stk_charge", "simu_energy_w_corr");
    auto h_stk_selected_charge2D = _fr_stk_charge_analysis.Filter("evtfilter_stk_charge_cut==true")
                                       .Histo2D<double, double, double>({"h_stk_selected_charge2D", "STK Charge", 500, 0, 100, 500, 0, 100}, "STK_chargeX", "STK_chargeY", "simu_energy_w_corr");

    // Extract NUD hisos
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_NUD_adc(DAMPE_NUD_channels);

    for (int channel = 0; channel < DAMPE_NUD_channels; ++channel)
        h_NUD_adc[channel] = _fr_bgo_analysis.Define("nud_adc_channel", [channel](std::vector<int> nud_adc) { return nud_adc[channel]; }, {"nud_adc"})
                                 .Histo1D<int, double>({(std::string("h_NUD_adc_") + std::to_string(channel)).c_str(), (std::string("NUD ADC - channel ") + std::to_string(channel)).c_str(), 100, 0, 1000}, "nud_adc_channel", "simu_energy_w_corr");
    auto h_NUD_total_adc = _fr_bgo_analysis.Histo1D<int, double>({"h_NUD_total_adc", "NUD Total ADC", 100, 0, 10000}, "nud_total_adc", "simu_energy_w_corr");
    auto h_NUD_max_adc = _fr_bgo_analysis.Histo1D<int, double>({"h_NUD_max_adc", "NUD Max ADC", 100, 0, 1000}, "nud_max_adc", "simu_energy_w_corr");
    auto h_NUD_max_channel = _fr_bgo_analysis.Histo1D<int, double>({"h_NUD_max_channel", "NUD Max Channel", 3, 0, 3}, "nud_max_channel_id", "simu_energy_w_corr");

    // Extract filter histos

    ROOT::RDF::RResultPtr<TH1D> h_geo_before_trigger_cut;
    if (_mc)
        h_geo_before_trigger_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                                        .Filter("evtfilter_geometric_before_trigger==true")
                                        .Histo1D<double, double>({"h_geo_before_trigger_cut", "Geometric factor pre trigger", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_trigger_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                         .Filter("evtfilter_evt_triggered==true")
                         .Histo1D<double, double>({"h_trigger_cut", "Trigger", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_geometric_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                           .Filter("evtfilter_geometric==true")
                           .Histo1D<double, double>({"h_geometric_cut", "Geometric", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_maxElayer_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                              .Filter("evtfilter_BGO_fiducial_maxElayer_cut==true")
                              .Histo1D<double, double>({"h_maxElayer_cut", "maxElayer cut", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_maxBarlayer_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                              .Filter("evtfilter_BGO_fiducial_maxBarLayer_cut==true")
                              .Histo1D<double, double>({"h_maxBarlayer_cut", "maxBarLayer cut", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_BGOTrackContainment_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                              .Filter("evtfilter_BGO_fiducial_BGOTrackContainment_cut==true")
                              .Histo1D<double, double>({"h_BGOTrackContainment_cut", "maxBarLayer cut", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_bgo_fiducial_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                              .Filter("evtfilter_BGO_fiducial==true")
                              .Histo1D<double, double>({"h_bgo_fiducial_cut", "BGO fiducial", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_nbarlayer13_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                              .Filter("evtfilter_nBarLayer13_cut==true")
                              .Histo1D<double, double>({"h_nbarlayer13_cut", "nBar layer 13 cut", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_maxrms_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                              .Filter("evtfilter_maxRms_cut==true")
                              .Histo1D<double, double>({"h_maxrms_cut", "max RMS cut", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_track_selection_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                                .Filter("evtfilter_track_selection_cut==true")
                                .Histo1D<double, double>({"h_track_selection_cut", "Track selection cut", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_psd_stk_match_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                                .Filter("evtfilter_psd_stk_match_cut==true")
                                .Histo1D<double, double>({"h_psd_stk_match_cut", "PSD-STK match cut", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_psd_charge_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                                .Filter("evtfilter_psd_charge_cut==true")
                                .Histo1D<double, double>({"h_psd_charge_cut", "PSD charge cut", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_stk_charge_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                                .Filter("evtfilter_stk_charge_cut==true")
                                .Histo1D<double, double>({"h_stk_charge_cut", "STK charge cut", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");
    auto h_all_cuts_cut = _fr_energy_filter.Define("energy_corr_gev", "energy_corr * 0.001")
                                .Filter("evtfilter_all_cut==true")
                                .Histo1D<double, double>({"h_all_cuts_cut", "All cuts", energy_nbins, &energy_binning[0]}, "energy_corr_gev", "simu_energy_w_corr");

    // Preselected events based histos

    auto h_BGOrec_ps_raw_energy = _fr_preselected.Define("raw_energy_gev", "energy * 0.001")
                                      .Histo1D<double, double>({"h_BGOrec_ps_raw_energy", "BGO Raw energy", energy_nbins, &energy_binning[0]}, "raw_energy_gev", "simu_energy_w_corr");
    auto h_BGOrec_ps_corr_energy = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                       .Histo1D<double, double>({"h_BGOrec_ps_corr_energy", "BGO Corr energy", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "simu_energy_w_corr");
    auto h_BGOrec_ps_cosine = _fr_preselected.Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                  .Histo1D<double, double>({"h_BGOrec_ps_cosine", "BGO Reco cosine - Preselection", 100, 0, 1}, "bgorec_cosine", "simu_energy_w_corr");

    auto h_BGOrec_ps_sumRms_bin_cosine2D_20_100 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                      .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                                      .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                      .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_bin_cosine2D_20_100", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_bin_cosine2D_100_250 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                       .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                                       .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                       .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_bin_cosine2D_100_250", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_bin_cosine2D_250_500 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                       .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                                       .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                       .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_bin_cosine2D_250_500", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_bin_cosine2D_500_1000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                        .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                                        .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                        .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_bin_cosine2D_500_1000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_bin_cosine2D_1000_3000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                         .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                                         .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                         .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_bin_cosine2D_1000_3000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_bin_cosine2D_3000_10000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                          .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                                          .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                          .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_bin_cosine2D_3000_10000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_bin_cosine2D_10000_20000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                           .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                                           .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                           .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_bin_cosine2D_10000_20000", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");

    auto h_BGOrec_ps_sumRms_flast_20_100 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                               .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                               .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_20_100", "sumRms vs F_{last} correlation - 20 GeV - 100 GeV; sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_100_250 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                                .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_100_250", "sumRms vs F_{last} correlation - 100 GeV - 250 GeV;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_250_500 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                                .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_250_500", "sumRms vs F_{last} correlation - 250 GeV - 500 GeV;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_500_1000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                 .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                                 .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_500_1000", "sumRms vs F_{last} correlation - 500 GeV - 1 TeV;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_1000_3000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                  .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                                  .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_1000_3000", "sumRms vs F_{last} correlation - 1 TeV - 3 TeV;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_3000_10000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                   .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                                   .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_3000_10000", "sumRms vs F_{last} correlation - 3 TeV - 10 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_10000_20000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                    .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                                    .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_10000_20000", "sumRms vs F_{last} correlation - 10 TeV - 20 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");

    auto h_BGOrec_ps_sumRms_flast_13_20_100 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                  .Filter("corr_energy_gev >= 20 && corr_energy_gev<100")
                                                  .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_13_20_100", "sumRms vs F_{13} correlation - 20 GeV - 100 GeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_13_100_250 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                   .Filter("corr_energy_gev >= 100 && corr_energy_gev<250")
                                                   .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_13_100_250", "sumRms vs F_{13} correlation - 100 GeV - 250 GeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_13_250_500 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                   .Filter("corr_energy_gev >= 250 && corr_energy_gev<500")
                                                   .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_13_250_500", "sumRms vs F_{13} correlation - 250 GeV - 500 GeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_13_500_1000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                    .Filter("corr_energy_gev >= 500 && corr_energy_gev<1000")
                                                    .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_13_500_1000", "sumRms vs F_{13} correlation - 500 GeV - 1 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_13_1000_3000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                     .Filter("corr_energy_gev >= 1000 && corr_energy_gev<3000")
                                                     .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_13_1000_3000", "sumRms vs F_{13} correlation - 1 TeV - 3 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_13_3000_10000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                      .Filter("corr_energy_gev >= 3000 && corr_energy_gev<10000")
                                                      .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_13_3000_10000", "sumRms vs F_{13} correlation - 3 TeV - 10 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_ps_sumRms_flast_13_10000_20000 = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                                                       .Filter("corr_energy_gev >= 10000 && corr_energy_gev<20000")
                                                       .Histo2D<double, double, double>({"h_BGOrec_ps_sumRms_flast_13_10000_20000", "sumRms vs F_{13} correlation - 10 TeV - 20 TeV ;sumRms [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
    auto h_BGOrec_ps_last_layer = _fr_preselected.Histo1D<int, double>({"h_BGOrec_ps_last_layer", "Last BGO layer; Last Energy Layer; entries", 14, 0, 14}, "lastBGOLayer", "simu_energy_w_corr");

    auto h_BGOrec_ps_ratio_last_layer = _fr_preselected.Histo1D<double, double>({"h_BGOrec_ps_ratio_last_layer", "Energy Ratio Last Layer", 200, flast_bin_lvalue, flast_bin_rvalue}, fracLast_leaf, "simu_energy_w_corr");

    std::shared_ptr<TH1D> h_BGOrec_ps_bar_energy = std::make_shared<TH1D>("h_BGOrec_ps_bar_energy", "BGO Bar Energy; Bar Energy [MeV]", 1e+3, 0, 1e+4);
    std::shared_ptr<TH1D> h_BGOrec_ps_bar_energy_zoom = std::make_shared<TH1D>("h_BGOrec_ps_bar_energy_zoom", "BGO Bar Energy; Bar Energy [MeV]", 1e+3, 0, 1e+2);
    std::shared_ptr<TH1D> h_BGOrec_ps_bar_energy_10MeV = std::make_shared<TH1D>("h_BGOrec_ps_bar_energy_10MeV", "BGO Bar Energy - 10 MeV filter; Bar Energy [MeV]", 1e+3, 10, 1e+4);
    std::shared_ptr<TH1D> h_BGOrec_ps_bar_energy_20MeV = std::make_shared<TH1D>("h_BGOrec_ps_bar_energy_20MeV", "BGO Bar Energy - 20 MeV filter; Bar Energy [MeV]", 1e+3, 20, 1e+4);
    std::shared_ptr<TH1D> h_BGOrec_ps_bar_energy_30MeV = std::make_shared<TH1D>("h_BGOrec_ps_bar_energy_30MeV", "BGO Bar Energy - 30 MeV filter; Bar Energy [MeV]", 1e+3, 30, 1e+4);
    std::shared_ptr<TH1D> h_BGOrec_ps_bar_energy_50MeV = std::make_shared<TH1D>("h_BGOrec_ps_bar_energy_50MeV", "BGO Bar Energy - 50 MeV filter; Bar Energy [MeV]", 1e+3, 50, 1e+4);
    std::shared_ptr<TH1D> h_BGOrec_ps_bar_energy_100MeV = std::make_shared<TH1D>("h_BGOrec_ps_bar_energy_100MeV", "BGO Bar Energy - 100 MeV filter; Bar Energy [MeV]", 1e+3, 100, 1e+4);

    auto h_BGOrec_ps_hits = _fr_preselected.Histo1D<int, double>({"h_BGOrec_ps_hits", "BGO Hits", 100, 0, 400}, "nBGOentries", "simu_energy_w_corr");
    std::shared_ptr<TH1D> h_BGOrec_ps_hits_10MeV = std::make_shared<TH1D>("h_BGOrec_ps_hits_10MeV", "BGO Hits - 10 MeV filter", 1e+2, 0, 4e+2);
    std::shared_ptr<TH1D> h_BGOrec_ps_hits_20MeV = std::make_shared<TH1D>("h_BGOrec_ps_hits_20MeV", "BGO Hits - 20 MeV filter", 1e+2, 0, 4e+2);
    std::shared_ptr<TH1D> h_BGOrec_ps_hits_30MeV = std::make_shared<TH1D>("h_BGOrec_ps_hits_30MeV", "BGO Hits - 30 MeV filter", 1e+2, 0, 4e+2);
    std::shared_ptr<TH1D> h_BGOrec_ps_hits_50MeV = std::make_shared<TH1D>("h_BGOrec_ps_hits_50MeV", "BGO Hits - 50 MeV filter", 1e+2, 0, 4e+2);
    std::shared_ptr<TH1D> h_BGOrec_ps_hits_100MeV = std::make_shared<TH1D>("h_BGOrec_ps_hits_100MeV", "BGO Hits - 100 MeV filter", 1e+2, 0, 4e+2);

    _fr_preselected.Foreach([
        &h_BGOrec_ps_bar_energy,
        &h_BGOrec_ps_bar_energy_zoom,
        &h_BGOrec_ps_bar_energy_10MeV, 
        &h_BGOrec_ps_bar_energy_20MeV, 
        &h_BGOrec_ps_bar_energy_30MeV, 
        &h_BGOrec_ps_bar_energy_50MeV, 
        &h_BGOrec_ps_bar_energy_100MeV,
        &h_BGOrec_ps_hits_10MeV,
        &h_BGOrec_ps_hits_20MeV,
        &h_BGOrec_ps_hits_30MeV,
        &h_BGOrec_ps_hits_50MeV,
        &h_BGOrec_ps_hits_100MeV]
        (std::vector<std::vector<double>> bar_energy, double energy_w) {
            unsigned int hits_10=0;
            unsigned int hits_20=0;
            unsigned int hits_30=0;
            unsigned int hits_50=0;
            unsigned int hits_100=0;

            for (auto& bgo_ly : bar_energy)
                for(auto& energy : bgo_ly)
                { 
                    h_BGOrec_ps_bar_energy->Fill(energy, energy_w);
                    if (energy <= 100)
                        h_BGOrec_ps_bar_energy_zoom->Fill(energy, energy_w);
                    if (energy>=10)
                    {
                        h_BGOrec_ps_bar_energy_10MeV->Fill(energy, energy_w);
                        ++hits_10;
                    }
                    if (energy>=20)
                    {
                        h_BGOrec_ps_bar_energy_20MeV->Fill(energy, energy_w);
                        ++hits_20;
                    }
                    if (energy>=30)
                    {
                        h_BGOrec_ps_bar_energy_30MeV->Fill(energy, energy_w);
                        ++hits_30;
                    }
                    if (energy>=50)
                    {
                        h_BGOrec_ps_bar_energy_50MeV->Fill(energy, energy_w);
                        ++hits_50;
                    }
                    if (energy>=100)
                    {
                        h_BGOrec_ps_bar_energy_100MeV->Fill(energy, energy_w);
                        ++hits_100;
                    }
                }
            
            h_BGOrec_ps_hits_10MeV->Fill(hits_10, energy_w);
            h_BGOrec_ps_hits_20MeV->Fill(hits_20, energy_w);
            h_BGOrec_ps_hits_30MeV->Fill(hits_30, energy_w);
            h_BGOrec_ps_hits_50MeV->Fill(hits_50, energy_w);
            h_BGOrec_ps_hits_100MeV->Fill(hits_100, energy_w);

        }, {"layerBarEnergy", "simu_energy_w_corr"});

    auto h_BGOrec_ps_slopeX = _fr_preselected.Histo1D<double, double>({"h_BGOrec_ps_slopeX", "BGOrec Slope X", 200, -10, 10}, "BGOrec_slopeX", "simu_energy_w_corr");
    auto h_BGOrec_ps_slopeY = _fr_preselected.Histo1D<double, double>({"h_BGOrec_ps_slopeY", "BGOrec Slope Y", 200, -10, 10}, "BGOrec_slopeY", "simu_energy_w_corr");
    auto h_BGOrec_ps_interceptX = _fr_preselected.Histo1D<double, double>({"h_BGOrec_ps_interceptX", "BGOrec Intercept X", 100, -500, 500}, "BGOrec_interceptX", "simu_energy_w_corr");
    auto h_BGOrec_ps_interceptY = _fr_preselected.Histo1D<double, double>({"h_BGOrec_ps_interceptY", "BGOrec Intercept Y", 100, -500, 500}, "BGOrec_interceptY", "simu_energy_w_corr");

    auto h_xtrl_ps_energy_int = _fr_preselected.Histo1D<double>({"h_xtrl_ps_energy_int", "XTRL", 200, 0, 150}, "xtrl", "simu_energy_w_corr");
    auto h_xtrl_ps = _fr_preselected.Define("corr_energy_gev", "energy_corr * 0.001")
                         .Histo2D<double, double, double>({"h_xtrl_ps", "XTRL", energy_nbins, &energy_binning[0], (int)xtrl_bins.size() - 1, &xtrl_bins[0]}, "corr_energy_gev", "xtrl", "simu_energy_w_corr");

    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ps_cosine_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ps_sumRms_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ps_sumRms_mean_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ps_sumRms_bin_weight(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_sumRms_sumRms_weight(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ps_sumRms_bin_cosine(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_sumRms_bin_cosine_cosine2D(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_sumRms_bin_cosine_cosine2D_log(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_sumRms_bin_cosine2D(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_sumRms_bin_cosine2D_log(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_sumRms_flast_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_sumRms_flast_13_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ps_max_energy_ratio(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ps_ratio_last(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ps_ratio_13(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ps_ratio_last_cosine(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_ratio_last_cosine_cosine2D(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_ratio_last_cosine_cosine2D_log(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_ratio_last_cosine2D(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_ratio_last_cosine2D_log(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_ratio_last_cosine2D_fdr(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_ratio_last_cosine2D_fdr_log(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ps_last_layer_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_BGOrec_ps_hits_bin(energy_nbins);

    std::vector<std::shared_ptr<TH1D>> h_BGOrec_ps_hits_bin_10MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_ps_hits_bin_20MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_ps_hits_bin_30MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_ps_hits_bin_50MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_ps_hits_bin_100MeV(energy_nbins);

    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_xtrl_ps_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_stkclusters_ps_bin(energy_nbins);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_stkclusters_bgohits_ps_bin(energy_nbins);
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH2D>>> h_stkclusters_ps_bgohits_layer(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH2D>>(DAMPE_stk_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_ps_rms_layer(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_ps_energy_ratio_layer(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_ps_energy_ratio_1R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_ps_energy_ratio_2R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_ps_energy_ratio_3R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_BGOrec_ps_energy_ratio_5R(energy_nbins, std::vector<ROOT::RDF::RResultPtr<TH1D>>(DAMPE_bgo_nLayers));

    std::vector<std::shared_ptr<TH1D>> h_BGOrec_ps_bar_energy_bin(energy_nbins);

    std::vector<std::shared_ptr<TH1D>> h_BGOrec_ps_bar_energy_bin_10_MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_ps_bar_energy_bin_20_MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_ps_bar_energy_bin_30_MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_ps_bar_energy_bin_50_MeV(energy_nbins);
    std::vector<std::shared_ptr<TH1D>> h_BGOrec_ps_bar_energy_bin_100_MeV(energy_nbins);

    std::vector<std::shared_ptr<TH2D>> h_BGOrec_ps_shower_profile(energy_nbins);
    std::vector<std::shared_ptr<TH2D>> h_BGOrec_ps_shower_profile_upto_09(energy_nbins);
    std::vector<std::shared_ptr<TH2D>> h_BGOrec_ps_shower_profile_from_09(energy_nbins);
    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        h_BGOrec_ps_bar_energy_bin[bin_idx - 1] = std::make_shared<TH1D>((std::string("h_BGOrec_ps_bar_energy_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 100, 0, 10000);
        h_BGOrec_ps_shower_profile[bin_idx - 1] = std::make_shared<TH2D>((std::string("h_BGOrec_ps_shower_profile_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Shower Profile - bin ") + std::to_string(bin_idx)).c_str(), 13, 0, 13, 100, 0, 1);
        h_BGOrec_ps_shower_profile_upto_09[bin_idx - 1] = std::make_shared<TH2D>((std::string("h_BGOrec_ps_shower_profile_upto_09_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Shower Profile cos(#theta) < 0.9 - bin ") + std::to_string(bin_idx)).c_str(), 13, 0, 13, 100, 0, 1);
        h_BGOrec_ps_shower_profile_from_09[bin_idx - 1] = std::make_shared<TH2D>((std::string("h_BGOrec_ps_shower_profile_from_09_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Shower Profile cos(#theta) > 0.9 - bin ") + std::to_string(bin_idx)).c_str(), 13, 0, 13, 100, 0, 1);

        h_BGOrec_ps_hits_bin_10MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_ps_hits_10MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Hits - 10 MeV filter - bin ") + std::to_string(bin_idx)).c_str(), 1e+2, 0, 4e+2);
        h_BGOrec_ps_hits_bin_20MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_ps_hits_20MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Hits - 20 MeV filter - bin ") + std::to_string(bin_idx)).c_str(), 1e+2, 0, 4e+2);
        h_BGOrec_ps_hits_bin_30MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_ps_hits_30MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Hits - 30 MeV filter - bin ") + std::to_string(bin_idx)).c_str(), 1e+2, 0, 4e+2);
        h_BGOrec_ps_hits_bin_50MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_ps_hits_50MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Hits - 50 MeV filter - bin ") + std::to_string(bin_idx)).c_str(), 1e+2, 0, 4e+2);
        h_BGOrec_ps_hits_bin_100MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_ps_hits_100MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Hits - 100 MeV filter - bin ") + std::to_string(bin_idx)).c_str(), 1e+2, 0, 4e+2);

        h_BGOrec_ps_bar_energy_bin_10_MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_ps_bar_energy_10_MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - 10 MeV filter - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 1e+3, 10, 1e+4);
        h_BGOrec_ps_bar_energy_bin_20_MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_ps_bar_energy_20_MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - 20 MeV filter - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 1e+3, 20, 1e+4);
        h_BGOrec_ps_bar_energy_bin_30_MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_ps_bar_energy_30_MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - 30 MeV filter - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 1e+3, 30, 1e+4);
        h_BGOrec_ps_bar_energy_bin_50_MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_ps_bar_energy_50_MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - 50 MeV filter - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 1e+3, 50, 1e+4);
        h_BGOrec_ps_bar_energy_bin_100_MeV[bin_idx -1] = std::make_shared<TH1D>((std::string("h_BGOrec_ps_bar_energy_100_MeV_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Bar Energy - 100 MeV filter - bin ") + std::to_string(bin_idx) + std::string("; Bar Energy [MeV]")).c_str(), 1e+3, 100, 1e+4);
    }

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        auto bin_filter = [=](int energy_bin) -> bool { return energy_bin == bin_idx; };
        h_BGOrec_ps_cosine_bin[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                  .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                  .Histo1D<double, double>({(std::string("h_BGOrec_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO Reco cosine - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "bgorec_cosine", "simu_energy_w_corr");
        h_BGOrec_ps_sumRms_bin[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                  .Histo1D<double, double>({(std::string("h_BGOrec_ps_sumRms_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms - bin ") + std::to_string(bin_idx)).c_str(), 100, sumrms_bin_lvalue, sumrms_bin_rvalue}, sumRms_leaf.c_str(), "simu_energy_w_corr");
        h_BGOrec_ps_sumRms_mean_bin[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                       .Define("num_layers", [](int last_layer) -> int { return last_layer + 1; }, {"lastBGOLayer"})
                                                       .Define("sumRms_eff", "sumRms/(lastBGOLayer+1)")
                                                       .Histo1D<double, double>({(std::string("h_BGOrec_ps_sumRms_mean_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms - bin ") + std::to_string(bin_idx)).c_str(), 80, 0, 100}, "sumRms_eff", "simu_energy_w_corr");
        h_BGOrec_ps_sumRms_bin_weight[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                         .Define("sumRms_w", computeSumRmsWeight, {"eLayer", "rmsLayer", "energy"})
                                                         .Histo1D<double, double>({(std::string("h_BGOrec_ps_sumRms_weight_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms weighted - bin ") + std::to_string(bin_idx)).c_str(), 50, 0, 100}, "sumRms_w", "simu_energy_w_corr");
        h_BGOrec_ps_sumRms_sumRms_weight[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                            .Define("sumRms_w", computeSumRmsWeight, {"eLayer", "rmsLayer", "energy"})
                                                            .Define("sumRms_eff", "sumRms/(lastBGOLayer+1)")
                                                            .Histo2D<double, double, double>({(std::string("h_BGOrec_ps_sumRms_sumRms_weight_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms_{mean} vs sumRms_{weighted} - bin ") + std::to_string(bin_idx) + std::string("; sumRms_{weighted} [mm]; sumRms_{mean} [mm]")).c_str(), 50, 0, 100, 50, 0, 100}, "sumRms_w", "sumRms_eff", "simu_energy_w_corr");
        h_BGOrec_ps_sumRms_bin_cosine[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                         .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                         .Define("sumRms_cosine", "sumRms/bgorec_cosine")
                                                         .Histo1D<double, double>({(std::string("h_BGOrec_ps_sumRms_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; sumRms/cos(#theta); counts")).c_str(), 100, 0, 3000}, "sumRms_cosine", "simu_energy_w_corr");
        h_BGOrec_ps_sumRms_bin_cosine_cosine2D[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                  .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                  .Define("sumRms_cosine", "sumRms/bgorec_cosine")
                                                                  .Histo2D<double, double, double>({(std::string("h_BGOrec_ps_sumRms_cosine_cosine2D_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); sumRms/cos(#theta)")).c_str(), 100, 0, 1, 100, 0, 3000}, "bgorec_cosine", "sumRms_cosine", "simu_energy_w_corr");
        h_BGOrec_ps_sumRms_bin_cosine_cosine2D_log[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                      .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                      .Define("sumRms_cosine", "sumRms/bgorec_cosine")
                                                                      .Histo2D<double, double, double>({(std::string("h_BGOrec_ps_sumRms_cosine_cosine2D_log_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); sumRms/cos(#theta)")).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_cosine_bins.size() - 1, &sumRms_cosine_bins[0]}, "bgorec_cosine", "sumRms_cosine", "simu_energy_w_corr");
        h_BGOrec_ps_sumRms_bin_cosine2D[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                           .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                           .Histo2D<double, double>({(std::string("h_BGOrec_ps_sumRms_cosine2D_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); sumRms [mm]")).c_str(), 100, 0, 1, 100, sumrms2D_bin_lvalue, sumrms2D_bin_rvalue}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
        if (!regularize_vars)
            h_BGOrec_ps_sumRms_bin_cosine2D_log[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                   .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                   .Histo2D<double, double>({(std::string("h_BGOrec_ps_sumRms_cosine2D_log_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); sumRms [mm]")).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)sumRms_binning.size() - 1, &sumRms_binning[0]}, "bgorec_cosine", sumRms_leaf.c_str(), "simu_energy_w_corr");
        h_BGOrec_ps_sumRms_flast_bin[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                        .Histo2D<double, double, double>({(std::string("h_BGOrec_ps_sumRms_flast_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms vs F_{last} - bin ") + std::to_string(bin_idx) + std::string("; sumRms [mm]; F_{last}")).c_str(), (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]}, sumRms_leaf.c_str(), fracLast_leaf, "simu_energy_w_corr");
        h_BGOrec_ps_sumRms_flast_13_bin[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                           .Histo2D<double, double, double>({(std::string("h_BGOrec_ps_sumRms_flast_13_bin_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms vs F_{13} - bin ") + std::to_string(bin_idx) + std::string("; sumRms [mm]; F_{13}")).c_str(), (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast13_binning.size() - 1, &flast13_binning[0]}, sumRms_leaf.c_str(), "fracLast_13", "simu_energy_w_corr");
        h_BGOrec_ps_max_energy_ratio[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                        .Define("max_energy_ratio", GetMaxEnergyRatio, {"fracLayer"})
                                                        .Histo1D<double, double>({(std::string("h_BGOrec_ps_max_energy_ratio_bin_") + std::to_string(bin_idx)).c_str(), (std::string("max energy ratio - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "max_energy_ratio", "simu_energy_w_corr");
        h_BGOrec_ps_ratio_last[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                  .Histo1D<double, double>({(std::string("h_BGOrec_ps_ratio_last_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer - bin ") + std::to_string(bin_idx)).c_str(), 200, flast_bin_lvalue, flast_bin_rvalue}, fracLast_leaf, "simu_energy_w_corr");
        h_BGOrec_ps_ratio_13[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                .Define("energy_ratio_13", "fracLayer[13]")
                                                .Histo1D<double, double>({(std::string("h_BGOrec_ps_ratio_13_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 13th Layer - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 0.01}, "energy_ratio_13", "simu_energy_w_corr");
        h_BGOrec_ps_ratio_last_cosine[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                         .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                         .Define("frac_last_cosine", "fracLast/bgorec_cosine")
                                                         .Histo1D<double, double>({(std::string("h_BGOrec_ps_ratio_last_cosine_bin_") + std::to_string(bin_idx)).c_str(), (std::string("F_{last}/cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; F_{last}/cos(#theta); entries")).c_str(), 1000, 0, 0.3}, "frac_last_cosine", "simu_energy_w_corr");
        h_BGOrec_ps_ratio_last_cosine_cosine2D[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                  .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                  .Define("frac_last_cosine", "fracLast/bgorec_cosine")
                                                                  .Histo2D<double, double, double>({(std::string("h_BGOrec_ps_ratio_last_cosine_cosine2D_bin_") + std::to_string(bin_idx)).c_str(), (std::string("F_{last}/cos(#theta) vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}/cos(#theta)")).c_str(), 100, 0, 1, 1000, 0, 0.3}, "bgorec_cosine", "frac_last_cosine", "simu_energy_w_corr");
        h_BGOrec_ps_ratio_last_cosine_cosine2D_log[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                      .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                      .Define("frac_last_cosine", "fracLast/bgorec_cosine")
                                                                      .Histo2D<double, double, double>({(std::string("h_BGOrec_ps_ratio_last_cosine_cosine2D_log_bin_") + std::to_string(bin_idx)).c_str(), (std::string("F_{last}/cos(#theta) vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}/cos(#theta)")).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)flast_cosine_binning.size() - 1, &flast_cosine_binning[0]}, "bgorec_cosine", "frac_last_cosine", "simu_energy_w_corr");
        h_BGOrec_ps_ratio_last_cosine2D[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                           .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                           .Histo2D<double, double, double>({(std::string("h_BGOrec_ps_ratio_last_cosine2D_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}")).c_str(), 100, 0, 1, 100, flast_bin_lvalue, flast_bin_rvalue}, "bgorec_cosine", fracLast_leaf, "simu_energy_w_corr");
        if (!regularize_vars)
        {
            h_BGOrec_ps_ratio_last_cosine2D_log[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                   .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                   .Histo2D<double, double, double>({(std::string("h_BGOrec_ps_ratio_last_cosine2D_log_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}")).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)flast_zoom_binning.size() - 1, &flast_zoom_binning[0]}, "bgorec_cosine", fracLast_leaf, "simu_energy_w_corr");
            h_BGOrec_ps_ratio_last_cosine2D_fdr[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                   .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                   .Histo2D<double, double, double>({(std::string("h_BGOrec_ps_ratio_last_cosine2D_fdr_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}")).c_str(), 100, 0, 1, 1000, 0, 0.2}, "bgorec_cosine", fracLast_leaf, "simu_energy_w_corr");
            h_BGOrec_ps_ratio_last_cosine2D_fdr_log[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                       .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                                       .Histo2D<double, double, double>({(std::string("h_BGOrec_ps_ratio_last_cosine2D_fdr_log_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio Last Layer vs cos(#theta) - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); F_{last}")).c_str(), (int)cosine_bins.size() - 1, &cosine_bins[0], (int)flast_binning.size() - 1, &flast_binning[0]}, "bgorec_cosine", fracLast_leaf, "simu_energy_w_corr");
        }
        h_BGOrec_ps_last_layer_bin[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                      .Histo1D<int, double>({(std::string("h_BGOrec_ps_last_layer_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Last BGO layer - bin ") + std::to_string(bin_idx) + std::string("; Last Energy Layer; entries")).c_str(), 14, 0, 14}, "lastBGOLayer", "simu_energy_w_corr");
        h_BGOrec_ps_hits_bin[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                            .Histo1D<int, double>({(std::string("h_BGOrec_ps_hits_bin_") + std::to_string(bin_idx)).c_str(), (std::string("BGO hits - bin ") + std::to_string(bin_idx)).c_str(), 80, 0, 350}, "nBGOentries", "simu_energy_w_corr");
        h_xtrl_ps_bin[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                         .Histo1D<double, double>({(std::string("h_xtrl_ps_bin_") + std::to_string(bin_idx)).c_str(), (std::string("XTRL - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 150}, "xtrl", "simu_energy_w_corr");

        h_stkclusters_ps_bin[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                        .Define("stk_clusters", [&bin_idx](const std::vector<int> plane_clusters) -> int { return std::accumulate(plane_clusters.begin(), plane_clusters.end(), 0); }, {"STK_plane_clusters"})
                                                .Histo1D<int, double>({(std::string("h_stkclusters_ps_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Total STK clusters - bin ") + std::to_string(bin_idx) + std::string(" ; total STK clusters; entries")).c_str(), 100, 0, 1500}, "stk_clusters", "simu_energy_w_corr");
        h_stkclusters_bgohits_ps_bin[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                .Define("stk_clusters", [](const std::vector<int> plane_clusters) -> int { return std::accumulate(plane_clusters.begin(), plane_clusters.end(), 0); }, {"STK_plane_clusters"})
                                                        .Histo2D<int, int, double>({(std::string("h_stkclusters_bgohist_ps_bin_") + std::to_string(bin_idx)).c_str(), (std::string("STK clusters vs BGO hits - bin ") + std::to_string(bin_idx) + std::string(" ; STK clusters; BGO hits")).c_str(), 100, 0, 500, 80, 0, 350}, "stk_clusters", "nBGOentries", "simu_energy_w_corr");

        _fr_preselected.Filter(bin_filter, {"energy_bin"})
            .Foreach([
                &h_BGOrec_ps_bar_energy_bin,
                &h_BGOrec_ps_bar_energy_bin_10_MeV,
                &h_BGOrec_ps_bar_energy_bin_20_MeV,
                &h_BGOrec_ps_bar_energy_bin_30_MeV,
                &h_BGOrec_ps_bar_energy_bin_50_MeV,
                &h_BGOrec_ps_bar_energy_bin_100_MeV,
                &h_BGOrec_ps_hits_bin_10MeV,
                &h_BGOrec_ps_hits_bin_20MeV,
                &h_BGOrec_ps_hits_bin_30MeV,
                &h_BGOrec_ps_hits_bin_50MeV,
                &h_BGOrec_ps_hits_bin_100MeV,
                bin_idx]
                (std::vector<std::vector<double>> bar_energy, double energy_w) 
                {
                    unsigned int hits_10=0;
                    unsigned int hits_20=0;
                    unsigned int hits_30=0;
                    unsigned int hits_50=0;
                    unsigned int hits_100=0;

                    for (auto& bgo_ly : bar_energy)
                        for(auto& energy : bgo_ly)
                        {
                            h_BGOrec_ps_bar_energy_bin[bin_idx -1]->Fill(energy, energy_w); 

                            if (energy>=10)
                            {
                                h_BGOrec_ps_bar_energy_bin_10_MeV[bin_idx -1]->Fill(energy, energy_w);
                                ++hits_10;
                            }
                            if (energy>=20)
                            {
                                h_BGOrec_ps_bar_energy_bin_20_MeV[bin_idx -1]->Fill(energy, energy_w);
                                ++hits_20;
                            }
                            if (energy>=30)
                            {
                                h_BGOrec_ps_bar_energy_bin_30_MeV[bin_idx -1]->Fill(energy, energy_w);
                                ++hits_30;
                            }
                            if (energy>=50)
                            {
                                h_BGOrec_ps_bar_energy_bin_50_MeV[bin_idx -1]->Fill(energy, energy_w);
                                ++hits_50;
                            }
                            if (energy>=100)
                            {
                                h_BGOrec_ps_bar_energy_bin_100_MeV[bin_idx -1]->Fill(energy, energy_w);
                                ++hits_50;
                            } 
                        }

                    h_BGOrec_ps_hits_bin_10MeV[bin_idx -1]->Fill(hits_10, energy_w);
                    h_BGOrec_ps_hits_bin_20MeV[bin_idx -1]->Fill(hits_20, energy_w);
                    h_BGOrec_ps_hits_bin_30MeV[bin_idx -1]->Fill(hits_30, energy_w);
                    h_BGOrec_ps_hits_bin_50MeV[bin_idx -1]->Fill(hits_50, energy_w);
                    h_BGOrec_ps_hits_bin_100MeV[bin_idx -1]->Fill(hits_100, energy_w);

                }, {"layerBarEnergy", "simu_energy_w_corr"});

        _fr_preselected.Filter(bin_filter, {"energy_bin"})
            .Foreach([&h_BGOrec_ps_shower_profile, bin_idx](std::vector<double> &energy_frac, double energy_w) { 
                for (int lidx=0; lidx<DAMPE_bgo_nLayers; ++lidx)
                   h_BGOrec_ps_shower_profile[bin_idx - 1]->Fill(lidx, energy_frac[lidx], energy_w); }, {"fracLayer", "simu_energy_w_corr"});

        _fr_preselected.Filter(bin_filter, {"energy_bin"})
            .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
            .Filter("bgorec_cosine<0.9")
            .Foreach([&h_BGOrec_ps_shower_profile_upto_09, bin_idx](std::vector<double> &energy_frac, double energy_w) { 
                for (int lidx=0; lidx<DAMPE_bgo_nLayers; ++lidx)
                    h_BGOrec_ps_shower_profile_upto_09[bin_idx - 1]->Fill(lidx, energy_frac[lidx], energy_w); }, {"fracLayer", "simu_energy_w_corr"});

        _fr_preselected.Filter(bin_filter, {"energy_bin"})
            .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
            .Filter("bgorec_cosine>0.9")
            .Foreach([&h_BGOrec_ps_shower_profile_from_09, bin_idx](std::vector<double> &energy_frac, double energy_w) { 
                for (int lidx=0; lidx<DAMPE_bgo_nLayers; ++lidx)
                    h_BGOrec_ps_shower_profile_from_09[bin_idx - 1]->Fill(lidx, energy_frac[lidx], energy_w); }, {"fracLayer", "simu_energy_w_corr"});

        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto GetLayerComponent = [=](std::vector<double> &value_layer) -> double { return value_layer[ly]; };
            auto GetStkLayerComponent = [=](std::vector<int> &stk_clusters) -> int { return stk_clusters[ly]; };

            if (ly < DAMPE_stk_nLayers)
                h_stkclusters_ps_bgohits_layer[bin_idx -1][ly] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                .Define("stk_layer_clusters", GetStkLayerComponent, {"STK_plane_clusters"})
                                                                .Histo2D<int, int, double>({(std::string("h_stkclusters_bgohits_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("STK clusters vs BGO hits - bin ") + std::to_string(bin_idx) + std::string(" ; STK clusters; BGO hits")).c_str(), 100, 0, 500, 80, 0, 350}, "stk_layer_clusters", "nBGOentries", "simu_energy_w_corr");
            h_BGOrec_ps_rms_layer[bin_idx - 1][ly] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                         .Define("rms_layer", GetLayerComponent, {"rmsLayer"})
                                                         .Histo1D<double, double>({(std::string("h_BGOrec_ps_rms_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("RMS layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 500}, "rms_layer", "simu_energy_w_corr");
            h_BGOrec_ps_energy_ratio_layer[bin_idx - 1][ly] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                                  .Define("energy_ratio_layer", GetLayerComponent, {"fracLayer"})
                                                                  .Histo1D<double, double>({(std::string("h_BGOrec_ps_energy_ratio_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_layer", "simu_energy_w_corr");
            h_BGOrec_ps_energy_ratio_1R[bin_idx - 1][ly] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                               .Define("energy_1R_layer", GetLayerComponent, {"energy_1R_radius"})
                                                               .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                               .Define("energy_ratio_1R_layer", "energy_1R_layer/energy_layer")
                                                               .Histo1D<double, double>({(std::string("h_BGOrec_ps_energy_ratio_1R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 1MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_1R_layer", "simu_energy_w_corr");
            h_BGOrec_ps_energy_ratio_2R[bin_idx - 1][ly] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                               .Define("energy_2R_layer", GetLayerComponent, {"energy_2R_radius"})
                                                               .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                               .Define("energy_ratio_2R_layer", "energy_2R_layer/energy_layer")
                                                               .Histo1D<double, double>({(std::string("h_BGOrec_ps_energy_ratio_2R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 2MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_2R_layer", "simu_energy_w_corr");
            h_BGOrec_ps_energy_ratio_3R[bin_idx - 1][ly] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                               .Define("energy_3R_layer", GetLayerComponent, {"energy_3R_radius"})
                                                               .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                               .Define("energy_ratio_3R_layer", "energy_3R_layer/energy_layer")
                                                               .Histo1D<double, double>({(std::string("h_BGOrec_ps_energy_ratio_3R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 3MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_3R_layer", "simu_energy_w_corr");
            h_BGOrec_ps_energy_ratio_5R[bin_idx - 1][ly] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                               .Define("energy_5R_layer", GetLayerComponent, {"energy_5R_radius"})
                                                               .Define("energy_layer", GetLayerComponent, {"eLayer"})
                                                               .Define("energy_ratio_5R_layer", "energy_5R_layer/energy_layer")
                                                               .Histo1D<double, double>({(std::string("h_BGOrec_ps_energy_ratio_5R_layer_") + std::to_string(ly) + std::string("_bin_") + std::to_string(bin_idx)).c_str(), (std::string("Energy Ratio 5MR layer ") + std::to_string(ly) + std::string(" - bin ") + std::to_string(bin_idx)).c_str(), 100, 0, 1}, "energy_ratio_5R_layer", "simu_energy_w_corr");
        }
    }

    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_energy;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_energy_w;
    ROOT::RDF::RResultPtr<TH1D> h_energy_ps_diff;
    ROOT::RDF::RResultPtr<TH1D> h_energy_ps_diff_corr;
    ROOT::RDF::RResultPtr<TH2D> h_energy_ps_diff2D;
    ROOT::RDF::RResultPtr<TH2D> h_energy_ps_diff2D_corr;
    ROOT::RDF::RResultPtr<TH2D> h_energy_ps_unfold;
    ROOT::RDF::RResultPtr<TH2D> h_energy_ps_unfold_corr;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_slopeX;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_slopeY;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_interceptX;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_interceptY;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_position_x;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_position_y;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_position_z;
    ROOT::RDF::RResultPtr<TH3D> h_simu_ps_position;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_flux_w;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_w;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_particle;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_theta;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_phi;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_charge;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_cosx;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_cosy;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_cosz;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_zenith;
    ROOT::RDF::RResultPtr<TH1D> h_simu_ps_azimuth;

    std::unique_ptr<TH1D> h_simu_ps_thruthtrajectory_start_x;
    std::unique_ptr<TH1D> h_simu_ps_thruthtrajectory_start_y;
    std::unique_ptr<TH1D> h_simu_ps_thruthtrajectory_start_z;
    std::unique_ptr<TH1D> h_simu_ps_thruthtrajectory_stop_x;
    std::unique_ptr<TH1D> h_simu_ps_thruthtrajectory_stop_y;
    std::unique_ptr<TH1D> h_simu_ps_thruthtrajectory_stop_z;
    std::unique_ptr<TH1D> h_simu_ps_thruthtrajectory_trackID;
    std::unique_ptr<TH1D> h_simu_ps_thruthtrajectory_parentID;
    std::unique_ptr<TH1D> h_simu_ps_thruthtrajectory_pdgID;
    std::unique_ptr<TH1D> h_simu_ps_thruthtrajectory_stop_index;

    if (_mc)
    {
        // Extract Simu histos
        h_simu_ps_energy = _fr_preselected.Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Histo1D<double, double>({"h_simu_ps_energy", "Simu energy", energy_nbins, &energy_binning[0]}, "simu_energy_gev", "simu_energy_w_corr");
        h_simu_ps_energy_w = _fr_preselected.Histo1D<double>({"h_simu_ps_energy_w", "Simu energy weight", 100, 0, 1}, "simu_energy_w_corr");
        h_energy_ps_diff = _fr_preselected.Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Define("raw_energy_gev", "energy * 0.001")
                                    .Define("energy_diff", "simu_energy_gev - raw_energy_gev")
                                    .Histo1D<double, double>({"h_energy_ps_diff", "Simu vs Raw Reco BGO energy: Simu Energy - Raw Energy (GeV); counts", 100, 0, 100}, "energy_diff", "simu_energy_w_corr");
        h_energy_ps_diff_corr = _fr_preselected.Define("simu_energy_gev", "simu_energy * 0.001")
                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Define("energy_diff", "simu_energy_gev - corr_energy_gev")
                                        .Histo1D<double, double>({"h_energy_ps_diff_corr", "Simu vs Corrected Reco BGO energy: Simu Energy - Corrected Energy (GeV); counts", 100, -100, 100}, "energy_diff", "simu_energy_w_corr");
        h_energy_ps_diff2D = _fr_preselected.Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Define("energy_ratio", "(energy - simu_energy)/simu_energy")
                                    .Histo2D<double, double, double>({"h_energy_ps_diff2D", "Energy Ratio; Real Energy (GeV); (Raw - Simu)/Simu", energy_nbins, &energy_binning[0], (int)energy_ratio_bins.size() - 1, &(energy_ratio_bins[0])}, "simu_energy_gev", "energy_ratio", "simu_energy_w_corr");
        h_energy_ps_diff2D_corr = _fr_preselected.Define("simu_energy_gev", "simu_energy * 0.001")
                                        .Define("energy_ratio", "(energy_corr-simu_energy)/simu_energy")
                                        .Histo2D<double, double, double>({"h_energy_ps_diff2D_corr", "Energy Ratio; Real Energy (GeV); (Corr - Simu)/Simu", energy_nbins, &energy_binning[0], (int)energy_ratio_bins.size() - 1, &(energy_ratio_bins[0])}, "simu_energy_gev", "energy_ratio", "simu_energy_w_corr");
        h_energy_ps_unfold = _fr_preselected.Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Define("raw_energy_gev", "energy * 0.001")
                                    .Histo2D<double, double, double>({"h_energy_ps_unfold", "Energy Unfolding Matrix; Real Energy (GeV); Raw Energy (GeV)", energy_nbins, &energy_binning[0], energy_nbins, &energy_binning[0]}, "simu_energy_gev", "raw_energy_gev", "simu_energy_w_corr");
        h_energy_ps_unfold_corr = _fr_preselected.Define("simu_energy_gev", "simu_energy * 0.001")
                                        .Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Histo2D<double, double, double>({"h_energy_ps_unfold_corr", "Energy Unfolding Matrix; Real Energy (GeV); Corr Energy (GeV)", energy_nbins, &energy_binning[0], energy_nbins, &energy_binning[0]}, "simu_energy_gev", "corr_energy_gev", "simu_energy_w_corr");
        h_simu_ps_slopeX = _fr_preselected.Histo1D<double, double>({"h_simu_ps_slopeX", "Simu Slope X", 100, -10, 10}, "simu_slope_x", "simu_energy_w_corr");
        h_simu_ps_slopeY = _fr_preselected.Histo1D<double, double>({"h_simu_ps_slopeY", "Simu Slope Y", 100, -10, 10}, "simu_slope_y", "simu_energy_w_corr");
        h_simu_ps_interceptX = _fr_preselected.Histo1D<double, double>({"h_simu_ps_interceptX", "Simu Intercept X", 100, -500, 500}, "simu_intercept_x", "simu_energy_w_corr");
        h_simu_ps_interceptY = _fr_preselected.Histo1D<double, double>({"h_simu_ps_interceptY", "Simu Intercept Y", 100, -500, 500}, "simu_intercept_y", "simu_energy_w_corr");
        h_simu_ps_position_x = _fr_preselected.Define("simu_position_comp", "simu_position.fX")
                                        .Histo1D<double, double>({"h_simu_ps_position_x", "Simu Position X; X [mm]", 100, -1500, 1500}, "simu_position_comp", "simu_energy_w_corr");
        h_simu_ps_position_y = _fr_preselected.Define("simu_position_comp", "simu_position.fY")
                                        .Histo1D<double, double>({"h_simu_ps_position_y", "Simu Position Y; Y [mm]", 100, -1500, 1500}, "simu_position_comp", "simu_energy_w_corr");
        h_simu_ps_position_z = _fr_preselected.Define("simu_position_comp", "simu_position.fZ")
                                        .Histo1D<double, double>({"h_simu_ps_position_z", "Simu Position Z; Z [mm]", 100, -1500, 1500}, "simu_position_comp", "simu_energy_w_corr");
        h_simu_ps_position = _fr_preselected.Define("simu_position_comp_x", "simu_position.fX")
                                    .Define("simu_position_comp_y", "simu_position.fY")
                                    .Define("simu_position_comp_z", "simu_position.fZ")
                                    .Histo3D<double, double, double, double>({"h_simu_ps_position", "Simu Position; Y [mm]; X [mm]; Z [mm]", 100, -1500, 1500, 100, -1500, 1500, 100, -1500, 1500}, "simu_position_comp_y", "simu_position_comp_x", "simu_position_comp_z", "simu_energy_w_corr");
        h_simu_ps_flux_w = _fr_preselected.Histo1D<double, double>({"h_simu_ps_flux_w", "Simu F_{w}; F_{w}", 100, 0, 2}, "simu_flux_w", "simu_energy_w_corr");
        h_simu_ps_w = _fr_preselected.Histo1D<double, double>({"h_simu_ps_w", "Simu w; w", 100, 0, 2}, "simu_w", "simu_energy_w_corr");
        h_simu_ps_particle = _fr_preselected.Histo1D<int, double>({"h_simu_ps_particle", "Simu particle", 3, 0, 3}, "simu_n_particle", "simu_energy_w_corr");
        h_simu_ps_theta = _fr_preselected.Histo1D<double, double>({"h_simu_ps_theta", "Simu #theta; #theta [deg]", 100, 0, 180}, "simu_theta", "simu_energy_w_corr");
        h_simu_ps_phi = _fr_preselected.Histo1D<double, double>({"h_simu_ps_phi", "Simu #phi; #phi [deg]", 200, -180, 180}, "simu_phi", "simu_energy_w_corr");
        h_simu_ps_charge = _fr_preselected.Histo1D<double, double>({"h_simu_ps_charge", "Simu charge;", 100, -2, 2}, "simu_charge", "simu_energy_w_corr");
        h_simu_ps_cosx = _fr_preselected.Histo1D<double, double>({"h_simu_ps_cosx", "Simu #cos(x);", 100, -1, 1}, "simu_cos_x", "simu_energy_w_corr");
        h_simu_ps_cosy = _fr_preselected.Histo1D<double, double>({"h_simu_ps_cosy", "Simu #cos(y);", 100, -1, 1}, "simu_cos_y", "simu_energy_w_corr");
        h_simu_ps_cosz = _fr_preselected.Histo1D<double, double>({"h_simu_ps_cosz", "Simu #cos(z);", 100, -1, 1}, "simu_cos_z", "simu_energy_w_corr");
        h_simu_ps_zenith = _fr_preselected.Histo1D<double, double>({"h_simu_ps_zenith", "Simu zenith; Zenith [deg]", 100, 0, 90}, "simu_zenith", "simu_energy_w_corr");
        h_simu_ps_azimuth = _fr_preselected.Histo1D<double, double>({"h_simu_ps_azimuth", "Simu azimuth; Azimuth [deg]", 200, -180, 180}, "simu_azimuth", "simu_energy_w_corr");

        h_simu_ps_thruthtrajectory_start_x = std::make_unique<TH1D>("h_simu_ps_thruthtrajectory_start_x", "Truth Trajectory Start X; X [mm]", 100, -1500, 1500);
        h_simu_ps_thruthtrajectory_start_y = std::make_unique<TH1D>("h_simu_ps_thruthtrajectory_start_y", "Truth Trajectory Start Y; Y [mm]", 100, -1500, 1500);
        h_simu_ps_thruthtrajectory_start_z = std::make_unique<TH1D>("h_simu_ps_thruthtrajectory_start_z", "Truth Trajectory Start Z; Z [mm]", 100, -1500, 1500);
        h_simu_ps_thruthtrajectory_stop_x = std::make_unique<TH1D>("h_simu_ps_thruthtrajectory_stop_x", "Truth Trajectory Stop X; X [mm]", 100, -3000, 3000);
        h_simu_ps_thruthtrajectory_stop_y = std::make_unique<TH1D>("h_simu_ps_thruthtrajectory_stop_y", "Truth Trajectory Stop Y; Y [mm]", 100, -3000, 3000);
        h_simu_ps_thruthtrajectory_stop_z = std::make_unique<TH1D>("h_simu_ps_thruthtrajectory_stop_z", "Truth Trajectory Stop Z; Z [mm]", 100, -3000, 3000);
        h_simu_ps_thruthtrajectory_trackID = std::make_unique<TH1D>("h_simu_ps_thruthtrajectory_trackID", "Truth Trajectory Track ID", 100, 0, 200);
        h_simu_ps_thruthtrajectory_parentID = std::make_unique<TH1D>("h_simu_ps_thruthtrajectory_parentID", "Truth Trajectory Parent ID", 2, 0, 1);
        h_simu_ps_thruthtrajectory_pdgID = std::make_unique<TH1D>("h_simu_ps_thruthtrajectory_pdgID", "Truth Trajectory PDG ID", 200, -2000, 3000);
        h_simu_ps_thruthtrajectory_stop_index = std::make_unique<TH1D>("h_simu_ps_thruthtrajectory_stop_index", "Truth Trajectory Stop Index", 100, -10, 10);

        _fr_preselected.Foreach([&h_simu_ps_thruthtrajectory_start_x](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_ps_thruthtrajectory_start_x->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_start_x", "simu_energy_w_corr"});
        _fr_preselected.Foreach([&h_simu_ps_thruthtrajectory_start_y](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_ps_thruthtrajectory_start_y->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_start_y", "simu_energy_w_corr"});
        _fr_preselected.Foreach([&h_simu_ps_thruthtrajectory_start_z](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_ps_thruthtrajectory_start_z->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_start_z", "simu_energy_w_corr"});
        _fr_preselected.Foreach([&h_simu_ps_thruthtrajectory_stop_x](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_ps_thruthtrajectory_stop_x->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_stop_x", "simu_energy_w_corr"});
        _fr_preselected.Foreach([&h_simu_ps_thruthtrajectory_stop_y](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_ps_thruthtrajectory_stop_y->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_stop_y", "simu_energy_w_corr"});
        _fr_preselected.Foreach([&h_simu_ps_thruthtrajectory_stop_z](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_ps_thruthtrajectory_stop_z->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_stop_z", "simu_energy_w_corr"});
        _fr_preselected.Foreach([&h_simu_ps_thruthtrajectory_trackID](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_ps_thruthtrajectory_trackID->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_trackID", "simu_energy_w_corr"});
        _fr_preselected.Foreach([&h_simu_ps_thruthtrajectory_parentID](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_ps_thruthtrajectory_parentID->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_parentID", "simu_energy_w_corr"});
        _fr_preselected.Foreach([&h_simu_ps_thruthtrajectory_pdgID](std::vector<double> thruthtrajectory, double energy_w) { 
                for (auto& _elm : thruthtrajectory)
                    h_simu_ps_thruthtrajectory_pdgID->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_PDG", "simu_energy_w_corr"});
        _fr_preselected.Foreach([&h_simu_ps_thruthtrajectory_stop_index](std::vector<double> thruthtrajectory, double energy_w) { 
            for (auto& _elm : thruthtrajectory)
                h_simu_ps_thruthtrajectory_stop_index->Fill(_elm, energy_w); }, {"simu_thruthtrajectory_stop_index", "simu_energy_w_corr"});
    }

    // Extract STK histos
    auto h_stk_ps_cosine = _fr_preselected.Histo1D({"h_stk_ps_cosine", "h_stk_cosine", 100, 0, 1}, "STK_bestTrack_costheta");
    auto h_stk_ps_slopeX = _fr_preselected.Histo1D<double, double>({"h_stk_ps_slopeX", "Simu Slope X", 200, -10, 10}, "STK_bestTrack_slopeX", "simu_energy_w_corr");
    auto h_stk_ps_slopeY = _fr_preselected.Histo1D<double, double>({"h_stk_ps_slopeY", "Simu Slope Y", 200, -10, 10}, "STK_bestTrack_slopeY", "simu_energy_w_corr");
    auto h_stk_ps_interceptX = _fr_preselected.Histo1D<double, double>({"h_stk_ps_interceptX", "Simu Intercept X", 100, -500, 500}, "STK_bestTrack_interceptX", "simu_energy_w_corr");
    auto h_stk_ps_interceptY = _fr_preselected.Histo1D<double, double>({"h_stk_ps_interceptY", "Simu Intercept Y", 100, -500, 500}, "STK_bestTrack_interceptY", "simu_energy_w_corr");

    // Extract PSD charge histos
    auto h_psd_ps_chargeX = _fr_preselected.Histo1D<double, double>({"h_psd_ps_chargeX", "PSD Charge X", 100, 0, 20}, "PSD_chargeX", "simu_energy_w_corr");
    auto h_psd_ps_chargeY = _fr_preselected.Histo1D<double, double>({"h_psd_ps_chargeY", "PSD Charge Y", 100, 0, 20}, "PSD_chargeY", "simu_energy_w_corr");
    auto h_psd_ps_charge = _fr_preselected.Define("psd_charge", [](double psd_chargeX, double psd_chargeY) { return 0.5 * (psd_chargeX + psd_chargeY); }, {"PSD_chargeX", "PSD_chargeY"})
                               .Histo1D<double, double>({"h_psd_ps_charge", "PSD Charge", 100, 0, 20}, "psd_charge", "simu_energy_w_corr");
    auto h_psd_ps_charge2D = _fr_preselected.Histo2D<double, double, double>({"h_psd_charge2D", "PSD Charge", 100, 0, 20, 100, 0, 20}, "PSD_chargeX", "PSD_chargeY", "simu_energy_w_corr");

    auto h_psd_ps_selected_chargeX = _fr_preselected.Filter("evtfilter_psd_charge_cut==true")
                                         .Histo1D<double, double>({"h_psd_ps_selected_chargeX", "PSD Charge X", 100, 0, 20}, "PSD_chargeX", "simu_energy_w_corr");
    auto h_psd_ps_selected_chargeY = _fr_preselected.Filter("evtfilter_psd_charge_cut==true")
                                         .Histo1D<double, double>({"h_psd_ps_selected_chargeY", "PSD Charge Y", 100, 0, 20}, "PSD_chargeY", "simu_energy_w_corr");
    auto h_psd_ps_selected_charge = _fr_preselected.Filter("evtfilter_psd_charge_cut==true")
                                        .Define("psd_charge", [](double psd_chargeX, double psd_chargeY) { return 0.5 * (psd_chargeX + psd_chargeY); }, {"PSD_chargeX", "PSD_chargeY"})
                                        .Histo1D<double, double>({"h_psd_ps_selected_charge", "PSD Charge", 100, 0, 20}, "psd_charge", "simu_energy_w_corr");
    auto h_psd_ps_selected_charge2D = _fr_preselected.Filter("evtfilter_psd_charge_cut==true")
                                          .Histo2D<double, double, double>({"h_psd_ps_selected_charge2D", "PSD Charge", 100, 0, 20, 100, 0, 20}, "PSD_chargeX", "PSD_chargeY", "simu_energy_w_corr");

    // Extract STK charge histos
    auto h_stk_ps_chargeX = _fr_preselected.Histo1D<double, double>({"h_stk_ps_chargeX", "STK Charge X", 100, 0, 20}, "STK_chargeX", "simu_energy_w_corr");
    auto h_stk_ps_chargeY = _fr_preselected.Histo1D<double, double>({"h_stk_ps_chargeY", "STK Charge Y", 100, 0, 20}, "STK_chargeY", "simu_energy_w_corr");
    auto h_stk_ps_charge = _fr_preselected.Define("stk_charge", [](double stk_chargeX, double stk_chargeY) { return 0.5 * (stk_chargeX + stk_chargeY); }, {"STK_chargeX", "STK_chargeY"})
                               .Histo1D<double, double>({"h_stk_charge", "STK Charge", 100, 0, 20}, "stk_charge", "simu_energy_w_corr");
    auto h_stk_ps_charge2D = _fr_preselected.Histo2D<double, double, double>({"h_stk_ps_charge2D", "STK Charge", 100, 0, 20, 100, 0, 20}, "STK_chargeX", "STK_chargeY", "simu_energy_w_corr");

    auto h_stk_ps_selected_chargeX = _fr_preselected.Filter("evtfilter_stk_charge_cut==true")
                                         .Histo1D<double, double>({"h_stk_ps_selected_chargeX", "STK Charge X", 100, 0, 20}, "STK_chargeX", "simu_energy_w_corr");
    auto h_stk_ps_selected_chargeY = _fr_preselected.Filter("evtfilter_stk_charge_cut==true")
                                         .Histo1D<double, double>({"h_stk_ps_selected_chargeY", "STK Charge Y", 100, 0, 20}, "STK_chargeY", "simu_energy_w_corr");
    auto h_stk_ps_selected_charge = _fr_preselected.Filter("evtfilter_stk_charge_cut==true")
                                        .Define("stk_charge", [](double stk_chargeX, double stk_chargeY) { return 0.5 * (stk_chargeX + stk_chargeY); }, {"STK_chargeX", "STK_chargeY"})
                                        .Histo1D<double, double>({"h_stk_ps_selected_charge", "STK Charge", 100, 0, 20}, "stk_charge", "simu_energy_w_corr");
    auto h_stk_ps_selected_charge2D = _fr_preselected.Filter("evtfilter_stk_charge_cut==true")
                                          .Histo2D<double, double, double>({"h_stk_ps_selected_charge2D", "STK Charge", 100, 0, 20, 100, 0, 20}, "STK_chargeX", "STK_chargeY", "simu_energy_w_corr");

    // Extract NUD hisos
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_NUD_ps_adc(DAMPE_NUD_channels);

    for (int channel = 0; channel < DAMPE_NUD_channels; ++channel)
        h_NUD_ps_adc[channel] = _fr_preselected.Define("nud_adc_channel", [channel](std::vector<int> nud_adc) { return nud_adc[channel]; }, {"nud_adc"})
                                    .Histo1D<int, double>({(std::string("h_NUD_ps_adc_") + std::to_string(channel)).c_str(), (std::string("NUD ADC - channel ") + std::to_string(channel)).c_str(), 100, 0, 1000}, "nud_adc_channel", "simu_energy_w_corr");
    auto h_NUD_ps_total_adc = _fr_preselected.Histo1D<int, double>({"h_NUD_ps_total_adc", "NUD Total ADC", 100, 0, 10000}, "nud_total_adc", "simu_energy_w_corr");
    auto h_NUD_ps_max_adc = _fr_preselected.Histo1D<int, double>({"h_NUD_ps_max_adc", "NUD Max ADC", 100, 0, 1000}, "nud_max_adc", "simu_energy_w_corr");
    auto h_NUD_ps_max_channel = _fr_preselected.Histo1D<int, double>({"h_NUD_ps_max_channel", "NUD Max Channel", 3, 0, 3}, "nud_max_channel_id", "simu_energy_w_corr");

    auto h_BGOrec_ps_sumRms_bin_cosine2D_complete = extractHistoFrorRegBH(energy_nbins);
    
    auto reg_sumrms_by_histo = [h_BGOrec_ps_sumRms_bin_cosine2D_complete] (const double sumrms, const double cosine, const int energy_bin) -> double {
        auto cosine_bin_idx = h_BGOrec_ps_sumRms_bin_cosine2D_complete[energy_bin-1]->GetXaxis()->FindBin(cosine);
        auto ptojy = static_cast<TH1D*>(h_BGOrec_ps_sumRms_bin_cosine2D_complete[energy_bin-1]->ProjectionY((std::string("h_BGOrec_ps_sumRms_bin_cosine2D_") + std::to_string(cosine_bin_idx)).c_str(), cosine_bin_idx, cosine_bin_idx));
        auto sumrms_mean_value = ptojy->GetMean();
        auto sumrms_rms_value = ptojy->GetRMS();
        auto sumrms_reg = sumrms_rms_value ? (sumrms-sumrms_mean_value)/sumrms_rms_value : sumrms-sumrms_mean_value;
        return sumrms_reg;
    };

    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_BGOrec_ps_sumRms_bin_cosine2D_hreg(energy_nbins);

    for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
    {
        auto bin_filter = [=](int energy_bin) -> bool { return energy_bin == bin_idx; };

        h_BGOrec_ps_sumRms_bin_cosine2D_hreg[bin_idx - 1] = _fr_preselected.Filter(bin_filter, {"energy_bin"})
                                                        .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                        .Define("bgorec_sumrms_reg", reg_sumrms_by_histo, {"sumRms", "bgorec_cosine", "energy_bin"})
                                                        .Histo2D<double, double>({(std::string("h_BGOrec_ps_sumRms_bin_cosine2D_hreg_") + std::to_string(bin_idx)).c_str(), (std::string("sumRms cosine - bin ") + std::to_string(bin_idx) + std::string("; cos(#theta); sumRms [mm]")).c_str(), 100, 0, 1, 100, -20, 20}, "bgorec_cosine", "bgorec_sumrms_reg", "simu_energy_w_corr");

    }

    if (verbose)
        std::cout << "Writing to disk... [" << outputPath << "]" << std::endl;

    TFile *outfile = TFile::Open(outputPath.c_str(), "RECREATE");

    if (_mc)
    {
        outfile->mkdir("Simu");
        outfile->cd("Simu");

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
        h_simu_position_x->Write();
        h_simu_position_y->Write();
        h_simu_position_z->Write();
        h_simu_position->Write();
        h_simu_flux_w->Write();
        h_simu_w->Write();
        h_simu_particle->Write();
        h_simu_theta->Write();
        h_simu_phi->Write();
        h_simu_charge->Write();
        h_simu_cosx->Write();
        h_simu_cosy->Write();
        h_simu_cosz->Write();
        h_simu_zenith->Write();
        h_simu_azimuth->Write();
        h_simu_thruthtrajectory_start_x->Write();
        h_simu_thruthtrajectory_start_y->Write();
        h_simu_thruthtrajectory_start_z->Write();
        h_simu_thruthtrajectory_stop_x->Write();
        h_simu_thruthtrajectory_stop_y->Write();
        h_simu_thruthtrajectory_stop_z->Write();
        h_simu_thruthtrajectory_trackID->Write();
        h_simu_thruthtrajectory_parentID->Write();
        h_simu_thruthtrajectory_pdgID->Write();
        h_simu_thruthtrajectory_stop_index->Write();
    }

    outfile->mkdir("STK");
    outfile->cd("STK");

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

    for (int bidx = 0; bidx < energy_nbins; ++bidx)
    {
        auto tmp_dir_name = std::string("STK/energybin_") + std::to_string(bidx + 1);
        outfile->mkdir(tmp_dir_name.c_str());
        outfile->cd(tmp_dir_name.c_str());
        h_stkclusters_bin[bidx]->Write();

        for (int ly = 0; ly < DAMPE_stk_nLayers; ++ly)
            h_stkclusters_bgohits_layer[bidx][ly]->Write();
    }

    outfile->mkdir("PSD");
    outfile->cd("PSD");

    h_psd_chargeX->Write();
    h_psd_chargeY->Write();
    h_psd_charge->Write();
    h_psd_charge2D->Write();
    h_psd_selected_chargeX->Write();
    h_psd_selected_chargeY->Write();
    h_psd_selected_charge->Write();
    h_psd_selected_charge2D->Write();

    outfile->mkdir("BGO");
    outfile->cd("BGO");

    h_BGOrec_raw_energy->Write();
    h_BGOrec_corr_energy->Write();
    h_BGOrec_cosine->Write();
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
    h_BGOrec_bar_energy->Write();
    h_BGOrec_bar_energy_zoom->Write();
    h_BGOrec_bar_energy_10MeV->Write();
    h_BGOrec_bar_energy_20MeV->Write();
    h_BGOrec_bar_energy_30MeV->Write();
    h_BGOrec_bar_energy_50MeV->Write();
    h_BGOrec_bar_energy_100MeV->Write();
    h_BGOrec_hits->Write();
    h_BGOrec_hits_10MeV->Write();
    h_BGOrec_hits_20MeV->Write();
    h_BGOrec_hits_30MeV->Write();
    h_BGOrec_hits_50MeV->Write();
    h_BGOrec_hits_100MeV->Write();
    h_BGOrec_last_layer->Write();

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
        h_BGOrec_sumRms_bin_cosine_cosine2D[bidx]->Write();
        h_BGOrec_sumRms_bin_cosine_cosine2D_log[bidx]->Write();
        h_BGOrec_sumRms_bin_cosine2D[bidx]->Write();

        if (!regularize_vars)
            h_BGOrec_sumRms_bin_cosine2D_log[bidx]->Write();
        h_BGOrec_sumRms_flast_bin[bidx]->Write();
        h_BGOrec_sumRms_flast_13_bin[bidx]->Write();
        h_BGOrec_max_energy_ratio[bidx]->Write();
        h_BGOrec_ratio_last[bidx]->Write();
        h_BGOrec_ratio_13[bidx]->Write();
        h_BGOrec_ratio_last_cosine[bidx]->Write();
        h_BGOrec_ratio_last_cosine_cosine2D[bidx]->Write();
        h_BGOrec_ratio_last_cosine_cosine2D_log[bidx]->Write();
        h_BGOrec_ratio_last_cosine2D[bidx]->Write();
        if (!regularize_vars)
        {
            h_BGOrec_ratio_last_cosine2D_log[bidx]->Write();
            h_BGOrec_ratio_last_cosine2D_fdr[bidx]->Write();
            h_BGOrec_ratio_last_cosine2D_fdr_log[bidx]->Write();
        }
        h_BGOrec_last_layer_bin[bidx]->Write();
        h_BGOrec_hits_bin[bidx]->Write();
        h_BGOrec_hits_bin_10MeV[bidx]->Write();
        h_BGOrec_hits_bin_20MeV[bidx]->Write();
        h_BGOrec_hits_bin_30MeV[bidx]->Write();
        h_BGOrec_hits_bin_50MeV[bidx]->Write();
        h_BGOrec_hits_bin_100MeV[bidx]->Write();
        h_xtrl_bin[bidx]->Write();
        h_stkclusters_bgohits_bin[bidx]->Write();
        h_BGOrec_bar_energy_bin[bidx]->Write();
        h_BGOrec_bar_energy_bin_10_MeV[bidx]->Write();
        h_BGOrec_bar_energy_bin_20_MeV[bidx]->Write();
        h_BGOrec_bar_energy_bin_30_MeV[bidx]->Write();
        h_BGOrec_bar_energy_bin_50_MeV[bidx]->Write();
        h_BGOrec_bar_energy_bin_100_MeV[bidx]->Write();
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

    outfile->mkdir("NUD");
    outfile->cd("NUD");

    for (auto &_elm : h_NUD_adc)
        _elm->Write();
    h_NUD_total_adc->Write();
    h_NUD_max_adc->Write();
    h_NUD_max_channel->Write();

    outfile->mkdir("Cuts");
    outfile->cd("Cuts");

    if (_mc)
        h_geo_before_trigger_cut->Write();
    h_trigger_cut->Write();
    h_geometric_cut->Write();
    h_maxElayer_cut->Write();
    h_maxBarlayer_cut->Write();
    h_BGOTrackContainment_cut->Write();
    h_bgo_fiducial_cut->Write();
    h_nbarlayer13_cut->Write();
    h_maxrms_cut->Write();
    h_track_selection_cut->Write();
    h_psd_stk_match_cut->Write();
    h_psd_charge_cut->Write();
    h_stk_charge_cut->Write();
    h_all_cuts_cut->Write();

    outfile->mkdir("Preselection");
    outfile->cd("Preselection");

    if (_mc)
    {
        outfile->mkdir("Preselection/Simu");
        outfile->cd("Preselection/Simu");

        h_simu_ps_energy->Write();
        h_simu_ps_energy_w->Write();
        h_energy_ps_diff->Write();
        h_energy_ps_diff_corr->Write();
        h_energy_ps_diff2D->Write();
        h_energy_ps_diff2D_corr->Write();
        h_energy_ps_unfold->Write();
        h_energy_ps_unfold_corr->Write();
        h_simu_ps_slopeX->Write();
        h_simu_ps_slopeY->Write();
        h_simu_ps_interceptX->Write();
        h_simu_ps_interceptY->Write();
        h_simu_ps_position_x->Write();
        h_simu_ps_position_y->Write();
        h_simu_ps_position_z->Write();
        h_simu_ps_position->Write();
        h_simu_ps_flux_w->Write();
        h_simu_ps_w->Write();
        h_simu_ps_particle->Write();
        h_simu_ps_theta->Write();
        h_simu_ps_phi->Write();
        h_simu_ps_charge->Write();
        h_simu_ps_cosx->Write();
        h_simu_ps_cosy->Write();
        h_simu_ps_cosz->Write();
        h_simu_ps_zenith->Write();
        h_simu_ps_azimuth->Write();
        h_simu_ps_thruthtrajectory_start_x->Write();
        h_simu_ps_thruthtrajectory_start_y->Write();
        h_simu_ps_thruthtrajectory_start_z->Write();
        h_simu_ps_thruthtrajectory_stop_x->Write();
        h_simu_ps_thruthtrajectory_stop_y->Write();
        h_simu_ps_thruthtrajectory_stop_z->Write();
        h_simu_ps_thruthtrajectory_trackID->Write();
        h_simu_ps_thruthtrajectory_parentID->Write();
        h_simu_ps_thruthtrajectory_pdgID->Write();
        h_simu_ps_thruthtrajectory_stop_index->Write();
    }

    outfile->mkdir("Preselection/STK");
    outfile->cd("Preselection/STK");

    h_stk_ps_cosine->Write();
    h_stk_ps_slopeX->Write();
    h_stk_ps_slopeY->Write();
    h_stk_ps_interceptX->Write();
    h_stk_ps_interceptY->Write();

    h_stk_ps_chargeX->Write();
    h_stk_ps_chargeY->Write();
    h_stk_ps_charge->Write();
    h_stk_ps_charge2D->Write();
    h_stk_ps_selected_chargeX->Write();
    h_stk_ps_selected_chargeY->Write();
    h_stk_ps_selected_charge->Write();
    h_stk_ps_selected_charge2D->Write();

    for (int bidx = 0; bidx < energy_nbins; ++bidx)
    {
        auto tmp_dir_name = std::string("Preselection/STK/energybin_") + std::to_string(bidx + 1);
        outfile->mkdir(tmp_dir_name.c_str());
        outfile->cd(tmp_dir_name.c_str());
        h_stkclusters_ps_bin[bidx]->Write();

        for (int ly = 0; ly < DAMPE_stk_nLayers; ++ly)
            h_stkclusters_ps_bgohits_layer[bidx][ly]->Write();
    }

    outfile->mkdir("Preselection/PSD");
    outfile->cd("Preselection/PSD");

    h_psd_ps_chargeX->Write();
    h_psd_ps_chargeY->Write();
    h_psd_ps_charge->Write();
    h_psd_ps_charge2D->Write();
    h_psd_ps_selected_chargeX->Write();
    h_psd_ps_selected_chargeY->Write();
    h_psd_ps_selected_charge->Write();
    h_psd_ps_selected_charge2D->Write();

    outfile->mkdir("Preselection/BGO");
    outfile->cd("Preselection/BGO");

    h_BGOrec_ps_raw_energy->Write();
    h_BGOrec_ps_corr_energy->Write();
    h_BGOrec_ps_cosine->Write();
    h_BGOrec_ps_slopeX->Write();
    h_BGOrec_ps_slopeY->Write();
    h_BGOrec_ps_interceptX->Write();
    h_BGOrec_ps_interceptY->Write();
    h_BGOrec_ps_sumRms_bin_cosine2D_20_100->Write();
    h_BGOrec_ps_sumRms_bin_cosine2D_100_250->Write();
    h_BGOrec_ps_sumRms_bin_cosine2D_250_500->Write();
    h_BGOrec_ps_sumRms_bin_cosine2D_500_1000->Write();
    h_BGOrec_ps_sumRms_bin_cosine2D_1000_3000->Write();
    h_BGOrec_ps_sumRms_bin_cosine2D_3000_10000->Write();
    h_BGOrec_ps_sumRms_bin_cosine2D_10000_20000->Write();
    h_BGOrec_ps_sumRms_flast_20_100->Write();
    h_BGOrec_ps_sumRms_flast_100_250->Write();
    h_BGOrec_ps_sumRms_flast_250_500->Write();
    h_BGOrec_ps_sumRms_flast_500_1000->Write();
    h_BGOrec_ps_sumRms_flast_1000_3000->Write();
    h_BGOrec_ps_sumRms_flast_3000_10000->Write();
    h_BGOrec_ps_sumRms_flast_10000_20000->Write();
    h_BGOrec_ps_sumRms_flast_13_20_100->Write();
    h_BGOrec_ps_sumRms_flast_13_100_250->Write();
    h_BGOrec_ps_sumRms_flast_13_250_500->Write();
    h_BGOrec_ps_sumRms_flast_13_500_1000->Write();
    h_BGOrec_ps_sumRms_flast_13_1000_3000->Write();
    h_BGOrec_ps_sumRms_flast_13_3000_10000->Write();
    h_BGOrec_ps_sumRms_flast_13_10000_20000->Write();
    h_xtrl_ps_energy_int->Write();
    h_xtrl_ps->Write();
    h_BGOrec_ps_bar_energy->Write();
    h_BGOrec_ps_bar_energy_zoom->Write();
    h_BGOrec_bar_energy_10MeV->Write();
    h_BGOrec_bar_energy_20MeV->Write();
    h_BGOrec_bar_energy_30MeV->Write();
    h_BGOrec_bar_energy_50MeV->Write();
    h_BGOrec_bar_energy_100MeV->Write();
    h_BGOrec_ps_hits->Write();
    h_BGOrec_ps_hits_10MeV->Write();
    h_BGOrec_ps_hits_20MeV->Write();
    h_BGOrec_ps_hits_30MeV->Write();
    h_BGOrec_ps_hits_50MeV->Write();
    h_BGOrec_ps_hits_100MeV->Write();
    h_BGOrec_ps_last_layer->Write();
    h_BGOrec_ps_ratio_last_layer->Write();

    for (int bidx = 0; bidx < energy_nbins; ++bidx)
    {
        auto tmp_dir_name = std::string("Preselection/BGO/energybin_") + std::to_string(bidx + 1);
        outfile->mkdir(tmp_dir_name.c_str());
        outfile->cd(tmp_dir_name.c_str());

        h_BGOrec_ps_cosine_bin[bidx]->Write();
        h_BGOrec_ps_sumRms_bin[bidx]->Write();
        h_BGOrec_ps_sumRms_mean_bin[bidx]->Write();
        h_BGOrec_ps_sumRms_bin_weight[bidx]->Write();
        h_BGOrec_ps_sumRms_sumRms_weight[bidx]->Write();
        h_BGOrec_ps_sumRms_bin_cosine[bidx]->Write();
        h_BGOrec_ps_sumRms_bin_cosine_cosine2D[bidx]->Write();
        h_BGOrec_ps_sumRms_bin_cosine_cosine2D_log[bidx]->Write();
        h_BGOrec_ps_sumRms_bin_cosine2D[bidx]->Write();
        h_BGOrec_ps_sumRms_bin_cosine2D_hreg[bidx]->Write();

        if (!regularize_vars)
            h_BGOrec_ps_sumRms_bin_cosine2D_log[bidx]->Write();
        h_BGOrec_ps_sumRms_flast_bin[bidx]->Write();
        h_BGOrec_ps_sumRms_flast_13_bin[bidx]->Write();
        h_BGOrec_ps_max_energy_ratio[bidx]->Write();
        h_BGOrec_ps_ratio_last[bidx]->Write();
        h_BGOrec_ps_ratio_13[bidx]->Write();
        h_BGOrec_ps_ratio_last_cosine[bidx]->Write();
        h_BGOrec_ps_ratio_last_cosine_cosine2D[bidx]->Write();
        h_BGOrec_ps_ratio_last_cosine_cosine2D_log[bidx]->Write();
        h_BGOrec_ps_ratio_last_cosine2D[bidx]->Write();
        if (!regularize_vars)
        {
            h_BGOrec_ps_ratio_last_cosine2D_log[bidx]->Write();
            h_BGOrec_ps_ratio_last_cosine2D_fdr[bidx]->Write();
            h_BGOrec_ps_ratio_last_cosine2D_fdr_log[bidx]->Write();
        }
        h_BGOrec_ps_last_layer_bin[bidx]->Write();
        h_BGOrec_ps_hits_bin[bidx]->Write();
        h_BGOrec_ps_hits_bin_10MeV[bidx]->Write();
        h_BGOrec_ps_hits_bin_20MeV[bidx]->Write();
        h_BGOrec_ps_hits_bin_30MeV[bidx]->Write();
        h_BGOrec_ps_hits_bin_50MeV[bidx]->Write();
        h_BGOrec_ps_hits_bin_100MeV[bidx]->Write();
        h_xtrl_ps_bin[bidx]->Write();
        h_stkclusters_bgohits_ps_bin[bidx]->Write();
        h_BGOrec_ps_bar_energy_bin[bidx]->Write();
        h_BGOrec_ps_bar_energy_bin_10_MeV[bidx]->Write();
        h_BGOrec_ps_bar_energy_bin_20_MeV[bidx]->Write();
        h_BGOrec_ps_bar_energy_bin_30_MeV[bidx]->Write();
        h_BGOrec_ps_bar_energy_bin_50_MeV[bidx]->Write();
        h_BGOrec_ps_bar_energy_bin_100_MeV[bidx]->Write();
        h_BGOrec_ps_shower_profile[bidx]->Write();
        h_BGOrec_ps_shower_profile_upto_09[bidx]->Write();
        h_BGOrec_ps_shower_profile_from_09[bidx]->Write();

        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            h_BGOrec_ps_rms_layer[bidx][ly]->Write();
            h_BGOrec_ps_energy_ratio_layer[bidx][ly]->Write();
            h_BGOrec_ps_energy_ratio_1R[bidx][ly]->Write();
            h_BGOrec_ps_energy_ratio_2R[bidx][ly]->Write();
            h_BGOrec_ps_energy_ratio_3R[bidx][ly]->Write();
            h_BGOrec_ps_energy_ratio_5R[bidx][ly]->Write();
        }
    }

    outfile->mkdir("Preselection/NUD");
    outfile->cd("Preselection/NUD");

    for (auto &_elm : h_NUD_ps_adc)
        _elm->Write();
    h_NUD_ps_total_adc->Write();
    h_NUD_ps_max_adc->Write();
    h_NUD_ps_max_channel->Write();
    
    outfile->Close();
}