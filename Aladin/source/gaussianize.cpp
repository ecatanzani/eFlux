#include "regularize.h"
#include "gaussianize.h"
#include "DAMPE_geo_structure.h"

#include <vector>
#include <cmath>
#include <map>

#include "TF1.h"
#include "TFile.h"
#include "TVector3.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <ROOT/RDataFrame.hxx>

void gaussianize(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<energy_config> _energy_config,
    std::shared_ptr<lambda_config> _lambda_config,
    const double _entries,
    const std::string outputPath,
    const std::string regularize_tree_path,
    const bool verbose,
    const unsigned int threads,
    const bool _mc)
{
    // Extract the energy info
    auto energy_binning = _energy_config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;
    auto min_evt_energy = _energy_config->GetMinEvtEnergy();
    auto max_evt_energy = _energy_config->GetMaxEvtEnergy();
    double _gev = 0.001;

    auto GetEnergyBin = [=](double energy) -> int  { 
        int bin_idx=0;
        for (; bin_idx<energy_nbins-1; ++bin_idx)
            if (energy * _gev >= energy_binning[bin_idx] && energy * _gev < energy_binning[bin_idx+1])
                break;
        return bin_idx+1; 
    };

    auto energyFilter = [&min_evt_energy, &max_evt_energy, &_gev](const double energy) -> bool {
        return energy*_gev >= min_evt_energy && energy*_gev <= max_evt_energy ? true : false;
    };

    auto simuEnergyWeight = [&_mc, &energy_binning](const double simu_energy_w) -> double {
        return _mc ? simu_energy_w * pow(energy_binning[0], 2) : 1;
    };

    // Create the RDF
    ROOT::EnableImplicitMT(threads);
    ROOT::RDataFrame _data_fr(*evtch);


    // Create the RDF with energy bin and correct energy weight
    auto _fr_energy_filter = _data_fr.Filter(energyFilter, {"energy_corr"})
                                        .Define("energy_bin", GetEnergyBin, {"energy_corr"})
                                        .Define("simu_energy_w_corr", simuEnergyWeight, {"simu_energy_w"});
    
    auto _fr_preselected = _fr_energy_filter.Filter("evtfilter_all_cut==true");

    std::cout << "\n\n**** Filter statistics ****\n";
    std::cout << "***************************\n";
    std::cout << "\nPreselected events: " << *(_fr_preselected.Count());
    std::cout << "\n\n********************";

    // Regularize SumRMS and ELF using angular distributions
    std::vector<TF1> sumrms_fitfunc(energy_nbins);
    std::vector<TF1> sumrms_fitfunc_err(energy_nbins);
    std::vector<TF1> flast_fitfunc(energy_nbins);
    std::vector<TF1> flast_fitfunc_err(energy_nbins);

    load_tf1s(
        regularize_tree_path,
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
    
    //*** Gaussianize variables

    auto rms_lambda_values = _lambda_config->GetRMSLambdaStruct();
    auto sumrms_lambda_values = _lambda_config->GetSumRMSLambdaStruct();
    auto elf_lambda_values = _lambda_config->GetELFLambdaStruct();
    auto ell_lambda_values = _lambda_config->GetELLLambdaStruct();
    auto xtrl_lambda_values = _lambda_config->GetXTRLLambdaStruct();

    if (verbose)
    {
        _lambda_config->PrintLambdaSettings();
        std::cout << "\nAnalysis running...\n";
    }

    auto gaussianize_rmslayer = [&rms_lambda_values](const std::vector<double> input_rmslayer) -> std::map<double, std::vector<double>>
    {
        auto gaussianize_elm = [](const std::vector<double> elm, const double lambda) -> std::vector<double> 
        {
            std::vector<double> elm_cp = elm;
            for (unsigned int idx=0; idx<elm.size(); ++idx)
                elm_cp[idx] = lambda ? (exp(lambda*elm[idx])-1)/lambda : elm[idx];
            return elm_cp;
        };
        
        std::map<double, std::vector<double>> rmslayer_gauss;
        for (int lambda_idx=0; lambda_idx<=rms_lambda_values.num; ++lambda_idx)
        {
            double lambda = rms_lambda_values.start + lambda_idx*rms_lambda_values.step;
            rmslayer_gauss.insert(std::pair<double, std::vector<double>>(lambda, gaussianize_elm(input_rmslayer, lambda)));
        }

        return rmslayer_gauss;
    };

    auto gaussianize_sumrms = [&sumrms_lambda_values](const double input_sumrmslayer) -> std::map<double, double>
    {
        auto gaussianize_elm = [](const double elm, const double lambda) -> double
        {
            double elm_cp = lambda ? (exp(lambda*elm)-1)/lambda : elm;
            return elm_cp;
        };
        
        std::map<double, double> sumrms_gauss;
        for (int lambda_idx=0; lambda_idx<=sumrms_lambda_values.num; ++lambda_idx)
        {
            double lambda = sumrms_lambda_values.start + lambda_idx*sumrms_lambda_values.step;
            sumrms_gauss.insert(std::pair<double, double>(lambda, gaussianize_elm(input_sumrmslayer, lambda)));
        }

        return sumrms_gauss;
    };

    auto gaussianize_fraclayer = [&elf_lambda_values](const std::vector<double> input_fraclayer) -> std::map<double, std::vector<double>>
    {
        auto gaussianize_elm = [](const std::vector<double> elm, const double lambda) -> std::vector<double> 
        {
            std::vector<double> elm_cp = elm;
            for (unsigned int idx=0; idx<elm.size(); ++idx)
                elm_cp[idx] = lambda ? (exp(lambda*elm[idx])-1)/lambda : elm[idx];
            return elm_cp;
        };

        std::map<double, std::vector<double>> fraclayer_gauss;
        for (int lambda_idx=0; lambda_idx<=elf_lambda_values.num; ++lambda_idx)
        {
            double lambda = elf_lambda_values.start + lambda_idx*elf_lambda_values.step;
            fraclayer_gauss.insert(std::pair<double, std::vector<double>>(lambda, gaussianize_elm(input_fraclayer, lambda)));
        }
        
        return fraclayer_gauss;
    };

    auto gaussianize_fraclastlayer = [&ell_lambda_values](const double input_fraclayer) -> std::map<double, double>
    {
        auto gaussianize_elm = [](const double elm, const double lambda) -> double
        {
            double elm_cp = lambda ? (exp(lambda*elm)-1)/lambda : elm;
            return elm_cp;
        };

        std::map<double, double> fraclastlayer_gauss;
        for (int lambda_idx=0; lambda_idx<=ell_lambda_values.num; ++lambda_idx)
        {
            double lambda = ell_lambda_values.start + lambda_idx*ell_lambda_values.step;
            fraclastlayer_gauss.insert(std::pair<double, double>(lambda, gaussianize_elm(input_fraclayer, lambda)));
        }
        
        return fraclastlayer_gauss;
    };

    auto gaussianize_xtrl = [&xtrl_lambda_values](const double input_xtrl) -> std::map<double, double>
    {
        auto gaussianize_elm = [](const double elm, const double lambda) -> double
        {
            double elm_cp = lambda ? (exp(lambda*elm)-1)/lambda : elm;
            return elm_cp;
        };

        std::map<double, double> xtrl_gauss;
        for (int lambda_idx=0; lambda_idx<=xtrl_lambda_values.num; ++lambda_idx)
        {
            double lambda = xtrl_lambda_values.start + lambda_idx*xtrl_lambda_values.step;
            xtrl_gauss.insert(std::pair<double, double>(lambda, gaussianize_elm(input_xtrl, lambda)));
        }
        
        return xtrl_gauss;
    };

    auto fr_gauss = fr.Define("rmslayer_gauss", gaussianize_rmslayer, {"rmsLayer"})
                        .Define("sumrms_gauss", gaussianize_sumrms, {"sumRms_reg"})
                        .Define("fraclayer_gauss", gaussianize_fraclayer, {"fracLayer"})
                        .Define("fraclastlayer_gauss", gaussianize_fraclastlayer, {"fracLast_reg"})
                        .Define("xtrl_gauss", gaussianize_xtrl, {"xtrl"});

    fr_gauss.Snapshot((std::string(evtch->GetName()) + std::string("_gauss")).c_str(), outputPath);
    
    if (verbose)
        std::cout << "\nOutput TFile has been written [" << outputPath << "]\n";
}