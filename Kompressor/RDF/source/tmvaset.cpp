#include "utils.h"
#include "tmvaset.h"
#include "binning.h"
#include "regularize.h"

#include "TF1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TVector3.h"
#include <ROOT/RDataFrame.hxx>

void createTMVAset(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    const std::string fit_tree_path,
    const bool signal,
    const std::string output_file,
    const bool _VERBOSE)
{
    // Enable multithreading
    ROOT::EnableImplicitMT();
    // Create RDF
    ROOT::RDataFrame _data_fr(*evtch);
    // Extract the energy binning
    auto energy_binning = _config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;
    // Check regularization option
    bool regularize_vars;
    double _gev = 0.001;
    fit_tree_path.empty() ? regularize_vars = false : regularize_vars = true;
    // Create the RDFs
    auto GetEnergyBin = [=](double energy) -> int { 
        int bin_idx=0;
        for (; bin_idx<energy_nbins-1; ++bin_idx)
            if (energy * _gev >= energy_binning[bin_idx] && energy * _gev < energy_binning[bin_idx+1])
                break;
        return bin_idx+1; };

    // Create the RDF with energy bin and Train/Test assignment variable
    auto _fr_bin_patch = _data_fr.Define("energy_bin", GetEnergyBin, {"energy_corr"})
                                .Define("tt_assign", []{return gRandom->Uniform();});

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
   
    if (_VERBOSE)
        std::cout << "\n\nINFO: Normalizing data facility has been activated\n";
    // Reading fit summary TTree
    load_tf1s(
        fit_tree_path,
        sumrms_fitfunc,
        sumrms_fitfunc_err,
        flast_fitfunc,
        flast_fitfunc_err,
        _VERBOSE);

    _fr_bin_patch = _fr_bin_patch.Define("sumRms_reg", regularize_sumrms, {"sumRms", "energy_bin", "BGOrec_trajectoryDirection2D"})
                        .Define("fracLast_reg", regularize_flast, {"fracLast", "energy_bin", "BGOrec_trajectoryDirection2D"});
    

    // Create the RDF for the preselected events
    auto _fr_preselected = _fr_bin_patch.Filter("evtfilter_all_cut==true");

    std::cout << "\n\n**** Filter statistics ****\n";
    std::cout << "***************************\n";
    std::cout << "\nPreselected events: " << *(_fr_preselected.Count());
    std::cout << "\n\n********************";

    // Create rms layer leafs
    _fr_preselected = _fr_preselected.Define("rmsLayer_1", "rmsLayer[0]")
                                    .Define("rmsLayer_2", "rmsLayer[1]")
                                    .Define("rmsLayer_3", "rmsLayer[2]")
                                    .Define("rmsLayer_4", "rmsLayer[3]")
                                    .Define("rmsLayer_5", "rmsLayer[4]")
                                    .Define("rmsLayer_6", "rmsLayer[5]")
                                    .Define("rmsLayer_7", "rmsLayer[6]")
                                    .Define("rmsLayer_8", "rmsLayer[7]")
                                    .Define("rmsLayer_9", "rmsLayer[8]")
                                    .Define("rmsLayer_10", "rmsLayer[9]")
                                    .Define("rmsLayer_11", "rmsLayer[10]")
                                    .Define("rmsLayer_12", "rmsLayer[11]")
                                    .Define("rmsLayer_13", "rmsLayer[12]")
                                    .Define("rmsLayer_14", "rmsLayer[13]");

    // Create flarLast layer leafs
    _fr_preselected = _fr_preselected.Define("fracLayer_1", "fracLayer[0]")
                                    .Define("fracLayer_2", "fracLayer[1]")
                                    .Define("fracLayer_3", "fracLayer[2]")
                                    .Define("fracLayer_4", "fracLayer[3]")
                                    .Define("fracLayer_5", "fracLayer[4]")
                                    .Define("fracLayer_6", "fracLayer[5]")
                                    .Define("fracLayer_7", "fracLayer[6]")
                                    .Define("fracLayer_8", "fracLayer[7]")
                                    .Define("fracLayer_9", "fracLayer[8]")
                                    .Define("fracLayer_10", "fracLayer[9]")
                                    .Define("fracLayer_11", "fracLayer[10]")
                                    .Define("fracLayer_12", "fracLayer[11]")
                                    .Define("fracLayer_13", "fracLayer[12]")
                                    .Define("fracLayer_14", "fracLayer[13]");
    
    // Create Train/Test frames
    auto _fr_train = _fr_preselected.Filter("tt_assign < 0.5");
    auto _fr_test = _fr_preselected.Filter("tt_assign > 0.5");
    
    // Write trees to file
    std::string train_tree_name;
    std::string test_tree_name;
    signal ? train_tree_name = "trainSignal" : train_tree_name = "trainBackground";
    signal ? test_tree_name = "testSignal" : test_tree_name = "testBackground";

    _fr_train.Snapshot(
        train_tree_name.c_str(), 
        expand_tt_output_path(output_file, 1).c_str(),
        {"STK_bestTrack_npoints", "energy", "energy_corr", "sumRms", "sumRms_reg", "fracLast", "fracLast_reg", 
        "rmsLayer_1", "rmsLayer_2", "rmsLayer_3", "rmsLayer_4", "rmsLayer_5", "rmsLayer_6", "rmsLayer_7", "rmsLayer_8", "rmsLayer_9", "rmsLayer_10", 
        "rmsLayer_11", "rmsLayer_12", "rmsLayer_13", "rmsLayer_14", "fracLayer_1", "fracLayer_2", "fracLayer_3", "fracLayer_4", "fracLayer_5","fracLayer_6", "fracLayer_7", "fracLayer_8", "fracLayer_9", "fracLayer_10", "fracLayer_11", "fracLayer_12", "fracLayer_13","fracLayer_14", "lastBGOLayer", "nBGOentries", "energy_1R_radius", "energy_2R_radius", "energy_3R_radius", "energy_5R_radius", "xtrl",
        "NUD_total_ADC.nud_total_adc", "NUD_max_ADC.nud_max_adc"});
    
    _fr_test.Snapshot(
        test_tree_name.c_str(), 
        expand_tt_output_path(output_file, 0).c_str(),
        {"STK_bestTrack_npoints", "energy", "energy_corr", "sumRms", "sumRms_reg", "fracLast", "fracLast_reg", 
        "rmsLayer_1", "rmsLayer_2", "rmsLayer_3", "rmsLayer_4", "rmsLayer_5", "rmsLayer_6", "rmsLayer_7", "rmsLayer_8", "rmsLayer_9", "rmsLayer_10", 
        "rmsLayer_11", "rmsLayer_12", "rmsLayer_13", "rmsLayer_14", "fracLayer_1", "fracLayer_2", "fracLayer_3", "fracLayer_4", "fracLayer_5","fracLayer_6", "fracLayer_7", "fracLayer_8", "fracLayer_9", "fracLayer_10", "fracLayer_11", "fracLayer_12", "fracLayer_13","fracLayer_14", "lastBGOLayer", "nBGOentries", "energy_1R_radius", "energy_2R_radius", "energy_3R_radius", "energy_5R_radius", "xtrl",
        "NUD_total_ADC.nud_total_adc", "NUD_max_ADC.nud_max_adc"});

}
