#include "tmva_vars.h"
#include "regularize.h"
#include "DAMPE_geo_structure.h"

#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TVector3.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <ROOT/RDataFrame.hxx>

struct best_lambda
{
    std::vector<std::vector<double>> rms;
    std::vector<double> sumrms;
    std::vector<std::vector<double>> fraclayer;
    std::vector<double> fraclast;
    std::vector<double> xtrl;

    std::vector<std::vector<double>> rms_norm_mean;
    std::vector<std::vector<double>> rms_norm_rms;
    std::vector<double> sumrms_norm_mean;
    std::vector<double> sumrms_norm_rms;
    std::vector<std::vector<double>> fraclayer_norm_mean;
    std::vector<std::vector<double>> fraclayer_norm_rms;
    std::vector<double> fraclast_norm_mean;
    std::vector<double> fraclast_norm_rms;
    std::vector<double> xtrl_norm_mean;
    std::vector<double> xtrl_norm_rms;

    void initSize(const int size)
    {
        rms.resize(size);
        sumrms.resize(size);
        fraclayer.resize(size);
        fraclast.resize(size);
        xtrl.resize(size);

        rms_norm_mean.resize(size);
        rms_norm_rms.resize(size);
        sumrms_norm_mean.resize(size);
        sumrms_norm_rms.resize(size);
        fraclayer_norm_mean.resize(size);
        fraclayer_norm_rms.resize(size);
        fraclast_norm_mean.resize(size);
        fraclast_norm_rms.resize(size);
        xtrl_norm_mean.resize(size);
        xtrl_norm_rms.resize(size);

        for (int idx=0; idx<size; ++idx)
        {
            rms[idx] = std::vector<double> (DAMPE_bgo_nLayers, 0);
            sumrms[idx] = 0;
            fraclayer[idx] = std::vector<double> (DAMPE_bgo_nLayers, 0);
            fraclast[idx] = 0;
            xtrl[idx] = 0;

            rms_norm_mean[idx] = std::vector<double> (DAMPE_bgo_nLayers, 0);
            rms_norm_rms[idx] = std::vector<double> (DAMPE_bgo_nLayers, 0);
            sumrms_norm_mean[idx] = 0;
            sumrms_norm_rms[idx] = 0;
            fraclayer_norm_mean[idx] = std::vector<double> (DAMPE_bgo_nLayers, 0);
            fraclayer_norm_rms[idx] = std::vector<double> (DAMPE_bgo_nLayers, 0);
            fraclast_norm_mean[idx] = 0;
            fraclast_norm_rms[idx] = 0;
            xtrl_norm_mean[idx] = 0;
            xtrl_norm_rms[idx] = 0;
        }

    };
};

inline best_lambda getLambdaBestValues(
    const std::string lambda_tree, 
    const int energy_nbins, 
    const bool verbose)
{
    TFile *lambda_tree_file = TFile::Open(lambda_tree.c_str(), "READ");
    if (!lambda_tree_file->IsOpen())
    {
        std::cerr << "\n\nError reading best lambda TTree [" << lambda_tree << "]\n\n";
        exit(100);
    }
    else
    {
        if (verbose)
            std::cout << "\n\nReading best lambda TTree [" << lambda_tree << "]";
    }

    auto lambda_select = best_lambda();
    lambda_select.initSize(energy_nbins);

    TTreeReader myReader("corrections_tree", lambda_tree_file);

    TTreeReaderValue<unsigned int> energy_bin(myReader, "energy_bin");
    TTreeReaderValue<std::vector<double>> best_rms_lambda(myReader, "best_rms_lambda");
    TTreeReaderValue<double> best_sumrms_lambda(myReader, "best_sumrms_lambda");
    TTreeReaderValue<std::vector<double>> best_fraclayer_lambda(myReader, "best_fraclayer_lambda");
    TTreeReaderValue<double> best_fraclast_lambda(myReader, "best_fraclast_lambda");
    TTreeReaderValue<double> best_xtrl_lambda(myReader, "best_xtrl_lambda");

    TTreeReaderValue<std::vector<double>> rms_norm_mean(myReader, "rms_norm_mean");
    TTreeReaderValue<std::vector<double>> rms_norm_rms(myReader, "rms_norm_rms");
    TTreeReaderValue<double> sumrms_norm_mean(myReader, "sumrms_norm_mean");
    TTreeReaderValue<double> sumrms_norm_rms(myReader, "sumrms_norm_rms");
    TTreeReaderValue<std::vector<double>> fraclayer_norm_mean(myReader, "fraclayer_norm_mean");
    TTreeReaderValue<std::vector<double>> fraclayer_norm_rms(myReader, "fraclayer_norm_rms");
    TTreeReaderValue<double> fraclast_norm_mean(myReader, "fraclast_norm_mean");
    TTreeReaderValue<double> fraclast_norm_rms(myReader, "fraclast_norm_rms");
    TTreeReaderValue<double> xtrl_norm_mean(myReader, "xtrl_norm_mean");
    TTreeReaderValue<double> xtrl_norm_rms(myReader, "xtrl_norm_rms");
    
    while (myReader.Next())
    {
        lambda_select.rms[*(energy_bin)-1] = *(best_rms_lambda);
        lambda_select.sumrms[*(energy_bin)-1] = *(best_sumrms_lambda);
        lambda_select.fraclayer[*(energy_bin)-1] = *(best_fraclayer_lambda);
        lambda_select.fraclast[*(energy_bin)-1] = *(best_fraclast_lambda);
        lambda_select.xtrl[*(energy_bin)-1] = *(best_xtrl_lambda);

        lambda_select.rms_norm_mean[*(energy_bin)-1] = *(rms_norm_mean);
        lambda_select.rms_norm_rms[*(energy_bin)-1] = *(rms_norm_rms);
        lambda_select.sumrms_norm_mean[*(energy_bin)-1] = *(sumrms_norm_mean);
        lambda_select.sumrms_norm_rms[*(energy_bin)-1] = *(sumrms_norm_rms);
        lambda_select.fraclayer_norm_mean[*(energy_bin)-1] = *(fraclayer_norm_mean);
        lambda_select.fraclast_norm_rms[*(energy_bin)-1] = *(fraclast_norm_rms);
        lambda_select.fraclast_norm_mean[*(energy_bin)-1] = *(fraclast_norm_mean);
        lambda_select.fraclast_norm_rms[*(energy_bin)-1] = *(fraclast_norm_rms);
        lambda_select.xtrl_norm_mean[*(energy_bin)-1] = *(xtrl_norm_mean);
        lambda_select.xtrl_norm_rms[*(energy_bin)-1] = *(xtrl_norm_rms);
    }

    lambda_tree_file->Close();
    
    return lambda_select;
}

void tmva_vars(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<energy_config> _energy_config,
    const double _entries,
    const std::string lambda_tree,
    const std::string regularize_tree,
    const std::string outputPath,
    const bool verbose,
    const unsigned int threads,
    const bool _mc)
{
    // Extract the energy binning
    auto energy_binning = _energy_config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;
    auto min_evt_energy = _energy_config->GetMinEvtEnergy();
    auto max_evt_energy = _energy_config->GetMaxEvtEnergy();
    double _gev = 0.001;

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

    // Create the RDF
    ROOT::EnableImplicitMT(threads);
    ROOT::RDataFrame _data_fr(*evtch);

    // Create the RDF with energy bin and correct energy weight
    auto _fr_energy_filter = _data_fr.Filter(energyFilter, {"energy_corr"})
                                        .Define("energy_bin", GetEnergyBin, {"energy_corr"});
    
    _fr_energy_filter = _mc ? _fr_energy_filter.Define("simu_energy_w_corr", simuEnergyWeight, {"simu_energy_w"}) : _fr_energy_filter.Define("simu_energy_w_corr", "1.").Define("simu_energy_w", "1.");

    auto _fr_preselected = _fr_energy_filter.Filter("evtfilter_all_cut==true");

    std::cout << "\n\n**** Filter statistics ****\n";
    std::cout << "***************************\n";
    std::cout << "\nPreselected events: " << *(_fr_preselected.Count());
    std::cout << "\n\n********************";

    auto lambda_select = getLambdaBestValues(lambda_tree, energy_nbins, verbose);

    // Regularize SumRMS and ELF using angular distributions
    std::vector<TF1> sumrms_fitfunc(energy_nbins);
    std::vector<TF1> sumrms_fitfunc_err(energy_nbins);
    std::vector<TF1> flast_fitfunc(energy_nbins);
    std::vector<TF1> flast_fitfunc_err(energy_nbins);

    load_tf1s(
        regularize_tree,
        sumrms_fitfunc,
        sumrms_fitfunc_err,
        flast_fitfunc,
        flast_fitfunc_err,
        verbose);

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
    
    auto fr = _fr_preselected.Define("sumRms_reg", regularize_sumrms, {"sumRms", "energy_bin", "BGOrec_trajectoryDirection2D"})
                            .Define("fracLast_reg", regularize_flast, {"fracLast", "energy_bin", "BGOrec_trajectoryDirection2D"});

    auto gaussianize_rmslayer = [lambda_select](const std::vector<double> input_rmslayer, const int energy_bin) -> std::vector<double>
    {
        auto gaussianize_elm = [](const std::vector<double> elm, const std::vector<double> lambda) -> std::vector<double> 
        {
            std::vector<double> elm_cp = elm;
            for (unsigned int idx=0; idx<elm.size(); ++idx)
                elm_cp[idx] = lambda[idx] ? (exp(lambda[idx]*elm[idx])-1)/lambda[idx] : elm[idx];
            return elm_cp;
        };
        
        return gaussianize_elm(input_rmslayer, lambda_select.rms[energy_bin-1]);
    };

    auto gaussianize_sumrms = [lambda_select](const double input_sumrmslayer, const int energy_bin) -> double
    {
        auto gaussianize_elm = [](const double elm, const double lambda) -> double
        {
            double elm_cp = lambda ? (exp(lambda*elm)-1)/lambda : elm;
            return elm_cp;
        };
        
        return gaussianize_elm(input_sumrmslayer, lambda_select.sumrms[energy_bin-1]);
    };

    auto gaussianize_fraclayer = [lambda_select](const std::vector<double> input_fraclayer, const int energy_bin) -> std::vector<double>
    {
        auto gaussianize_elm = [](const std::vector<double> elm, const std::vector<double> lambda) -> std::vector<double> 
        {
            std::vector<double> elm_cp = elm;
            for (unsigned int idx=0; idx<elm.size(); ++idx)
                elm_cp[idx] = lambda[idx] ? (exp(lambda[idx]*elm[idx])-1)/lambda[idx] : elm[idx];
            return elm_cp;
        };
        
        return gaussianize_elm(input_fraclayer, lambda_select.fraclayer[energy_bin-1]);
    };

    auto gaussianize_fraclastlayer = [lambda_select](const double input_fraclayer, const int energy_bin) -> double
    {
        auto gaussianize_elm = [](const double elm, const double lambda) -> double
        {
            double elm_cp = lambda ? (exp(lambda*elm)-1)/lambda : elm;
            return elm_cp;
        };

        return gaussianize_elm(input_fraclayer, lambda_select.fraclast[energy_bin-1]);
    };
    
    auto gaussianize_xtrl = [lambda_select](const double input_xtrl, const int energy_bin) -> double
    {
        auto gaussianize_elm = [](const double elm, const double lambda) -> double
        {
            double elm_cp = lambda ? (exp(lambda*elm)-1)/lambda : elm;
            return elm_cp;
        };

        return gaussianize_elm(input_xtrl, lambda_select.xtrl[energy_bin-1]);
    };

    auto fr_gauss = fr.Define("rmslayer_gauss", gaussianize_rmslayer, {"rmsLayer", "energy_bin"})
                        .Define("sumrms_gauss", gaussianize_sumrms, {"sumRms_reg", "energy_bin"})
                        .Define("fraclayer_gauss", gaussianize_fraclayer, {"fracLayer", "energy_bin"})
                        .Define("fraclastlayer_gauss", gaussianize_fraclastlayer, {"fracLast_reg", "energy_bin"})
                        .Define("xtrl_gauss", gaussianize_xtrl, {"xtrl", "energy_bin"});

    auto normalize_rmslayer = [lambda_select](const std::vector<double> rms_gauss, const int energy_bin) -> std::vector<double>
    {
        std::vector<double> rms_norm = rms_gauss;
        for (int l_idx=0; l_idx<DAMPE_bgo_nLayers; ++l_idx)
        {
            if (lambda_select.rms_norm_rms[energy_bin-1][l_idx]) 
                rms_norm[l_idx] /= lambda_select.rms_norm_rms[energy_bin-1][l_idx];
            if (lambda_select.rms_norm_mean[energy_bin-1][l_idx]>0)
                rms_norm[l_idx] -= lambda_select.rms_norm_mean[energy_bin-1][l_idx];
            else
                rms_norm[l_idx] += abs(lambda_select.rms_norm_mean[energy_bin-1][l_idx]);
        }
        return rms_norm;
    };

    auto normalize_sumrms = [lambda_select](const double sumrms_gauss, const int energy_bin) -> double
    {
        double sumrms_norm = sumrms_gauss;
        if (lambda_select.sumrms_norm_rms[energy_bin-1])
            sumrms_norm /= lambda_select.sumrms_norm_rms[energy_bin-1];
        if (lambda_select.sumrms_norm_mean[energy_bin-1]>0)
            sumrms_norm -= lambda_select.sumrms_norm_rms[energy_bin-1];
        else
            sumrms_norm += abs(lambda_select.sumrms_norm_rms[energy_bin-1]);
        return sumrms_norm;
    };

    auto normalize_fraclayer = [lambda_select](const std::vector<double> fraclayer_gauss, const int energy_bin) -> std::vector<double>
    {
        std::vector<double> fraclayer_norm = fraclayer_gauss;
        for (int l_idx=0; l_idx<DAMPE_bgo_nLayers; ++l_idx)
        {
            if (lambda_select.fraclayer_norm_rms[energy_bin-1][l_idx]) 
                fraclayer_norm[l_idx] /= lambda_select.fraclayer_norm_rms[energy_bin-1][l_idx];
            if (lambda_select.fraclayer_norm_mean[energy_bin-1][l_idx]>0)
                fraclayer_norm[l_idx] -= lambda_select.fraclayer_norm_mean[energy_bin-1][l_idx];
            else
                fraclayer_norm[l_idx] += abs(lambda_select.fraclayer_norm_mean[energy_bin-1][l_idx]);
        }
        return fraclayer_norm;
    };

    auto normalize_fraclastlayer = [lambda_select](const double fraclastlayer_gauss, const int energy_bin) -> double
    {
        double fraclastlayer_norm = fraclastlayer_gauss;
        if (lambda_select.fraclast_norm_rms[energy_bin-1])
            fraclastlayer_norm /= lambda_select.fraclast_norm_rms[energy_bin-1];
        if (lambda_select.sumrms_norm_mean[energy_bin-1]>0)
            fraclastlayer_norm -= lambda_select.fraclast_norm_mean[energy_bin-1];
        else
            fraclastlayer_norm += abs(lambda_select.fraclast_norm_mean[energy_bin-1]);
        return fraclastlayer_norm;
    };

    auto normalize_xtrl = [lambda_select](const double xtrl_gauss, const int energy_bin) -> double
    {
        double xtrl_norm = xtrl_gauss;
        if (lambda_select.xtrl_norm_rms[energy_bin-1])
            xtrl_norm /= lambda_select.xtrl_norm_rms[energy_bin-1];
        if (lambda_select.xtrl_norm_mean[energy_bin-1]>0)
            xtrl_norm -= lambda_select.xtrl_norm_mean[energy_bin-1];
        else
            xtrl_norm += abs(lambda_select.xtrl_norm_mean[energy_bin-1]);
        return xtrl_norm;
    };

    auto fr_norm = fr_gauss.Define("rmslayer_norm", normalize_rmslayer, {"rmslayer_gauss", "energy_bin"})
                        .Define("sumrms_norm", normalize_sumrms, {"sumrms_gauss", "energy_bin"})
                        .Define("fraclayer_norm", normalize_fraclayer, {"fraclayer_gauss", "energy_bin"})
                        .Define("fraclastlayer_norm", normalize_fraclastlayer, {"fraclastlayer_gauss", "energy_bin"})
                        .Define("xtrl_norm", normalize_xtrl, {"xtrl_gauss", "energy_bin"});

    auto fr_tmva = fr_norm.Define("rmslayer_norm_1", "rmslayer_norm[0]")
                            .Define("rmslayer_norm_2", "rmslayer_norm[1]")
                            .Define("rmslayer_norm_3", "rmslayer_norm[2]")
                            .Define("rmslayer_norm_4", "rmslayer_norm[3]")
                            .Define("rmslayer_norm_5", "rmslayer_norm[4]")
                            .Define("rmslayer_norm_6", "rmslayer_norm[5]")
                            .Define("rmslayer_norm_7", "rmslayer_norm[6]")
                            .Define("rmslayer_norm_8", "rmslayer_norm[7]")
                            .Define("rmslayer_norm_9", "rmslayer_norm[8]")
                            .Define("rmslayer_norm_10", "rmslayer_norm[9]")
                            .Define("rmslayer_norm_11", "rmslayer_norm[10]")
                            .Define("rmslayer_norm_12", "rmslayer_norm[11]")
                            .Define("rmslayer_norm_13", "rmslayer_norm[12]")
                            .Define("rmslayer_norm_14", "rmslayer_norm[13]")

                            .Define("fraclayer_norm_1", "fraclayer_norm[0]")
                            .Define("fraclayer_norm_2", "fraclayer_norm[1]")
                            .Define("fraclayer_norm_3", "fraclayer_norm[2]")
                            .Define("fraclayer_norm_4", "fraclayer_norm[3]")
                            .Define("fraclayer_norm_5", "fraclayer_norm[4]")
                            .Define("fraclayer_norm_6", "fraclayer_norm[5]")
                            .Define("fraclayer_norm_7", "fraclayer_norm[6]")
                            .Define("fraclayer_norm_8", "fraclayer_norm[7]")
                            .Define("fraclayer_norm_9", "fraclayer_norm[8]")
                            .Define("fraclayer_norm_10", "fraclayer_norm[9]")
                            .Define("fraclayer_norm_11", "fraclayer_norm[10]")
                            .Define("fraclayer_norm_12", "fraclayer_norm[11]")
                            .Define("fraclayer_norm_13", "fraclayer_norm[12]")
                            .Define("fraclayer_norm_14", "fraclayer_norm[13]");

    fr_tmva.Snapshot((std::string(evtch->GetName()) + std::string("_norm")).c_str(), outputPath);
    
    if (verbose)
        std::cout << "\nOutput TFile has been written [" << outputPath << "]\n";
}