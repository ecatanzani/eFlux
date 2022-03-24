#include "vars.h"

#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"

inline double xtrl_computation(const double sumRms, const double lastFracLayer)
{
	return lastFracLayer != -1 ? 0.125e-6 * pow(sumRms, 4) * lastFracLayer : -999;
}

void vars::load_cosine_corrections(std::string cosine_regularize_path, const bool verbose)
{
    TFile *fitfile = TFile::Open(cosine_regularize_path.c_str(), "READ");
    if (!fitfile->IsOpen())
    {
        std::cerr << "\n\nError reading summary fit TTree [" << cosine_regularize_path << "]\n\n";
        exit(100);
    }
    else
        if (verbose)
            std::cout << "\nReading fitting summary [" << cosine_regularize_path << "]\n";

    // **** Fit function parameters
    const char* tree_name = "fit_summary";
    const char* func_form = "pol3";
    const int npars = 4;
    const double lvalue = 0;
    const double rvalue = 1; 
    // *****

    auto my_corrections_tree = static_cast<TTree*>(fitfile->Get(tree_name));
    std::vector<double>* flast_pars         {nullptr};
    std::vector<double>* flast_err_pars     {nullptr};
    std::vector<double>* sumrms_pars        {nullptr};
    std::vector<double>* sumrms_err_pars    {nullptr};
    my_corrections_tree->SetBranchAddress("flast_pars", &flast_pars);
    my_corrections_tree->SetBranchAddress("flast_err_pars", &flast_err_pars);
    my_corrections_tree->SetBranchAddress("sumrms_pars", &sumrms_pars);
    my_corrections_tree->SetBranchAddress("sumrms_err_pars", &sumrms_err_pars);

    energy_nbins = my_corrections_tree->GetEntries();
    cosine_corrections.sumrms_fitfunc.resize(energy_nbins);
    cosine_corrections.sumrms_fitfunc_err.resize(energy_nbins);
    cosine_corrections.flast_fitfunc.resize(energy_nbins);
    cosine_corrections.flast_fitfunc_err.resize(energy_nbins);

    for (int bin_idx=0; bin_idx<energy_nbins; ++bin_idx)
    {
        my_corrections_tree->GetEntry(bin_idx);

        cosine_corrections.sumrms_fitfunc[bin_idx]      = TF1(
            (std::string(tree_name) + std::string("_sumrms_fitfunc_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);
        cosine_corrections.sumrms_fitfunc_err[bin_idx]  = TF1(
            (std::string(tree_name) + std::string("_sumrms_fitfunc_err_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);
        cosine_corrections.flast_fitfunc[bin_idx]       = TF1(
            (std::string(tree_name) + std::string("_flast_fitfunc_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);
        cosine_corrections.flast_fitfunc_err[bin_idx]   = TF1(
            (std::string(tree_name) + std::string("_flast_fitfunc_err_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);

        for (int par_idx=0; par_idx<npars; ++par_idx)
        {
            cosine_corrections.sumrms_fitfunc[bin_idx]      .SetParameter(par_idx, sumrms_pars->at(par_idx));
            cosine_corrections.sumrms_fitfunc_err[bin_idx]  .SetParameter(par_idx, sumrms_err_pars->at(par_idx));
            cosine_corrections.flast_fitfunc[bin_idx]       .SetParameter(par_idx, flast_pars->at(par_idx));
            cosine_corrections.flast_fitfunc_err[bin_idx]   .SetParameter(par_idx, flast_err_pars->at(par_idx));
        }
        
    }

    fitfile->Close();
}

void vars::load_box_cox_corrections(std::string box_cox_regularize_path, const bool verbose)
{
    TFile *lambda_tree_file = TFile::Open(box_cox_regularize_path.c_str(), "READ");
    if (!lambda_tree_file->IsOpen())
    {
        std::cerr << "\n\nError reading best lambda TTree [" << box_cox_regularize_path << "]\n\n";
        exit(100);
    }
    else
    {
        if (verbose)
            std::cout << "\n\nReading best lambda TTree [" << box_cox_regularize_path << "]";
    }

    box_cox_correction_parameters.initSize(energy_nbins);
    auto my_corrections_tree = static_cast<TTree*>(lambda_tree_file->Get("corrections_tree"));
    
    unsigned int energy_bin;
    std::vector<double>* best_rms_lambda             {nullptr};
    std::vector<double>* best_fraclayer_lambda       {nullptr};
    std::vector<double>* rms_norm_mean               {nullptr};
    std::vector<double>* rms_norm_rms                {nullptr};
    std::vector<double>* fraclayer_norm_mean         {nullptr};
    std::vector<double>* fraclayer_norm_rms          {nullptr};
    
    double best_sumrms_lambda;
    double best_fraclast_lambda;
    double best_xtrl_lambda;
    double sumrms_norm_mean;
    double sumrms_norm_rms;
    double fraclast_norm_mean;
    double fraclast_norm_rms;
    double xtrl_norm_mean;
    double xtrl_norm_rms;

    my_corrections_tree->SetBranchAddress("energy_bin", &energy_bin);
    my_corrections_tree->SetBranchAddress("best_rms_lambda", &best_rms_lambda);
    my_corrections_tree->SetBranchAddress("best_fraclayer_lambda", &best_fraclayer_lambda);
    my_corrections_tree->SetBranchAddress("rms_norm_mean", &rms_norm_mean);
    my_corrections_tree->SetBranchAddress("rms_norm_rms", &rms_norm_rms);
    my_corrections_tree->SetBranchAddress("fraclayer_norm_mean", &fraclayer_norm_mean);
    my_corrections_tree->SetBranchAddress("fraclayer_norm_rms", &fraclayer_norm_rms);
    my_corrections_tree->SetBranchAddress("best_sumrms_lambda", &best_sumrms_lambda);
    my_corrections_tree->SetBranchAddress("best_fraclast_lambda", &best_fraclast_lambda);
    my_corrections_tree->SetBranchAddress("best_xtrl_lambda", &best_xtrl_lambda);
    my_corrections_tree->SetBranchAddress("sumrms_norm_mean", &sumrms_norm_mean);
    my_corrections_tree->SetBranchAddress("sumrms_norm_rms", &sumrms_norm_rms);
    my_corrections_tree->SetBranchAddress("fraclast_norm_mean", &fraclast_norm_mean);
    my_corrections_tree->SetBranchAddress("fraclast_norm_rms", &fraclast_norm_rms);
    my_corrections_tree->SetBranchAddress("xtrl_norm_mean", &xtrl_norm_mean);
    my_corrections_tree->SetBranchAddress("xtrl_norm_rms", &xtrl_norm_rms);

    for (int bin_idx=0; bin_idx<energy_nbins; ++bin_idx)
    {
        my_corrections_tree->GetEntry(bin_idx);

        box_cox_correction_parameters.rms[energy_bin-1]                     = *best_rms_lambda;
        box_cox_correction_parameters.sumrms[energy_bin-1]                  = best_sumrms_lambda;
        box_cox_correction_parameters.fraclayer[energy_bin-1]               = *best_fraclayer_lambda;
        box_cox_correction_parameters.fraclast[energy_bin-1]                = best_fraclast_lambda;
        box_cox_correction_parameters.xtrl[energy_bin-1]                    = best_xtrl_lambda;

        box_cox_correction_parameters.rms_norm_mean[energy_bin-1]           = *rms_norm_mean;
        box_cox_correction_parameters.rms_norm_rms[energy_bin-1]            = *rms_norm_rms;

        box_cox_correction_parameters.sumrms_norm_mean[energy_bin-1]        = sumrms_norm_mean;
        box_cox_correction_parameters.sumrms_norm_rms[energy_bin-1]         = sumrms_norm_rms;

        box_cox_correction_parameters.fraclayer_norm_mean[energy_bin-1]     = *fraclayer_norm_mean;
        box_cox_correction_parameters.fraclayer_norm_rms[energy_bin-1]      = *fraclayer_norm_rms;

        box_cox_correction_parameters.fraclast_norm_mean[energy_bin-1]      = fraclast_norm_mean;
        box_cox_correction_parameters.fraclast_norm_rms[energy_bin-1]       = fraclast_norm_rms;

        box_cox_correction_parameters.xtrl_norm_mean[energy_bin-1]          = xtrl_norm_mean;
        box_cox_correction_parameters.xtrl_norm_rms[energy_bin-1]           = xtrl_norm_rms;
    }

    lambda_tree_file->Close();
}

vars::vars(
    const std::string cosine_regularize_path,
	const std::string box_cox_regularize_path,
	const bool verbose)
    {
        // Load cosine angular corrections
        load_cosine_corrections(cosine_regularize_path, verbose);
        // Load box-cox lambda corrections
        load_box_cox_corrections(box_cox_regularize_path, verbose);
    }

const bdt_vars vars::GetVars(
    const std::vector<double> &rms,
    const double sumrms,
    const std::vector<double> &fraclayer,
    const double fraclastlayer,
    const double corrected_energy_gev,
    const std::vector<float> &energy_binning,
    const TVector3& bgo_direction)
    {
        // Reset the variables struct
        my_vars.Reset();
        
        // Initialize with the new event values
        double xtrl = xtrl_computation(sumrms, fraclastlayer);
        my_vars.corrected_energy_gev = corrected_energy_gev;
        my_vars.xtrl_spectator = xtrl;
        
        // Apply cosine regularization
        auto get_energy_bin = [&energy_binning] (const double corrected_energy) -> int
        {
            int energybin {1};
            for (size_t idx=0; idx<energy_binning.size()-1; ++idx)
                if (corrected_energy>=energy_binning[idx] && corrected_energy<energy_binning[idx+1])
                    energybin = idx+1;
            return energybin;
        };

        auto regularize_sumrms = [=](double sumrms, int energy_bin, TVector3 bgodir) -> double 
        {
            // Initialize regularized sumrms variable
            double reg_sumrms = sumrms;
            // Initialize BGO cosine from directrion
            double bgocosine = bgodir.CosTheta();
            // Regularize sumrms
            reg_sumrms -= cosine_corrections.sumrms_fitfunc[energy_bin - 1].Eval(bgocosine);
            reg_sumrms /= cosine_corrections.sumrms_fitfunc_err[energy_bin - 1].Eval(bgocosine);
            return reg_sumrms;
        };

        auto regularize_flast = [=](double flast, int energy_bin, TVector3 bgodir) -> double 
        {
            // Initialize regularized sumrms variable
            double reg_flast = flast;
            // Initialize BGO cosine from directrion
            double bgocosine = bgodir.CosTheta();
            // Regularize sumrms
            reg_flast -= cosine_corrections.flast_fitfunc[energy_bin - 1].Eval(bgocosine);
            reg_flast /= cosine_corrections.flast_fitfunc_err[energy_bin - 1].Eval(bgocosine);
            return reg_flast;
        };

        my_vars.sumrms = regularize_sumrms(sumrms, get_energy_bin(corrected_energy_gev), bgo_direction);
        my_vars.fraclastlayer = regularize_flast(fraclastlayer, get_energy_bin(corrected_energy_gev), bgo_direction);
        
        // Apply box-cox regularization
        auto gaussianize_rmslayer = [=](const std::vector<double> input_rmslayer, const int energy_bin) -> std::vector<double>
        {
            auto gaussianize_elm = [](const std::vector<double> elm, const std::vector<double> lambda) -> std::vector<double> 
            {
                std::vector<double> elm_cp (elm.size(), 0);
                for (unsigned int idx=0; idx<elm.size(); ++idx)
                    elm_cp[idx] = lambda[idx] ? (exp(lambda[idx]*elm[idx])-1)/lambda[idx] : elm[idx];
                return elm_cp;
            };
            
            return gaussianize_elm(input_rmslayer, box_cox_correction_parameters.rms[energy_bin-1]);
        };

        auto gaussianize_sumrms = [=](const double input_sumrmslayer, const int energy_bin) -> double
        {
            auto gaussianize_elm = [](const double elm, const double lambda) -> double
            {
                double elm_cp = lambda ? (exp(lambda*elm)-1)/lambda : elm;
                return elm_cp;
            };
            
            return gaussianize_elm(input_sumrmslayer, box_cox_correction_parameters.sumrms[energy_bin-1]);
        };

        auto gaussianize_fraclayer = [=](const std::vector<double> input_fraclayer, const int energy_bin) -> std::vector<double>
        {
            auto gaussianize_elm = [](const std::vector<double> elm, const std::vector<double> lambda) -> std::vector<double> 
            {
                std::vector<double> elm_cp (elm.size(), 0);
                for (unsigned int idx=0; idx<elm.size(); ++idx)
                    elm_cp[idx] = lambda[idx] ? (exp(lambda[idx]*elm[idx])-1)/lambda[idx] : elm[idx];
                return elm_cp;
            };
            
            return gaussianize_elm(input_fraclayer, box_cox_correction_parameters.fraclayer[energy_bin-1]);
        };

        auto gaussianize_fraclastlayer = [=](const double input_fraclayer, const int energy_bin) -> double
        {
            auto gaussianize_elm = [](const double elm, const double lambda) -> double
            {
                double elm_cp = lambda ? (exp(lambda*elm)-1)/lambda : elm;
                return elm_cp;
            };

            return gaussianize_elm(input_fraclayer, box_cox_correction_parameters.fraclast[energy_bin-1]);
        };
        
        auto gaussianize_xtrl = [=](const double input_xtrl, const int energy_bin) -> double
        {
            auto gaussianize_elm = [](const double elm, const double lambda) -> double
            {
                double elm_cp = lambda ? (exp(lambda*elm)-1)/lambda : elm;
                return elm_cp;
            };

            return gaussianize_elm(input_xtrl, box_cox_correction_parameters.xtrl[energy_bin-1]);
        };

        my_vars.rms = gaussianize_rmslayer(rms, get_energy_bin(corrected_energy_gev));
        my_vars.sumrms = gaussianize_sumrms(my_vars.sumrms, get_energy_bin(corrected_energy_gev));
        my_vars.fraclayer = gaussianize_fraclayer(fraclayer, get_energy_bin(corrected_energy_gev));
        my_vars.fraclastlayer = gaussianize_fraclastlayer(my_vars.fraclastlayer, get_energy_bin(corrected_energy_gev));
        my_vars.xtrl = gaussianize_xtrl(xtrl, get_energy_bin(corrected_energy_gev));
        
        auto normalize_rmslayer = [=](const std::vector<double> rms_gauss, const int energy_bin) -> std::vector<double>
        {
            std::vector<double> rms_norm = rms_gauss;
            for (int l_idx=0; l_idx<DAMPE_bgo_nLayers; ++l_idx)
            {
                if (box_cox_correction_parameters.rms_norm_rms[energy_bin-1][l_idx]) 
                    rms_norm[l_idx] /= box_cox_correction_parameters.rms_norm_rms[energy_bin-1][l_idx];
                if (box_cox_correction_parameters.rms_norm_mean[energy_bin-1][l_idx]>0)
                    rms_norm[l_idx] -= box_cox_correction_parameters.rms_norm_mean[energy_bin-1][l_idx];
                else
                    rms_norm[l_idx] += abs(box_cox_correction_parameters.rms_norm_mean[energy_bin-1][l_idx]);
            }
            return rms_norm;
        };

        auto normalize_sumrms = [=](const double sumrms_gauss, const int energy_bin) -> double
        {
            double sumrms_norm = sumrms_gauss;
            if (box_cox_correction_parameters.sumrms_norm_rms[energy_bin-1])
                sumrms_norm /= box_cox_correction_parameters.sumrms_norm_rms[energy_bin-1];
            if (box_cox_correction_parameters.sumrms_norm_mean[energy_bin-1]>0)
                sumrms_norm -= box_cox_correction_parameters.sumrms_norm_rms[energy_bin-1];
            else
                sumrms_norm += abs(box_cox_correction_parameters.sumrms_norm_rms[energy_bin-1]);
            return sumrms_norm;
        };

        auto normalize_fraclayer = [=](const std::vector<double> fraclayer_gauss, const int energy_bin) -> std::vector<double>
        {
            std::vector<double> fraclayer_norm = fraclayer_gauss;
            for (int l_idx=0; l_idx<DAMPE_bgo_nLayers; ++l_idx)
            {
                if (box_cox_correction_parameters.fraclayer_norm_rms[energy_bin-1][l_idx]) 
                    fraclayer_norm[l_idx] /= box_cox_correction_parameters.fraclayer_norm_rms[energy_bin-1][l_idx];
                if (box_cox_correction_parameters.fraclayer_norm_mean[energy_bin-1][l_idx]>0)
                    fraclayer_norm[l_idx] -= box_cox_correction_parameters.fraclayer_norm_mean[energy_bin-1][l_idx];
                else
                    fraclayer_norm[l_idx] += abs(box_cox_correction_parameters.fraclayer_norm_mean[energy_bin-1][l_idx]);
            }
            return fraclayer_norm;
        };

        auto normalize_fraclastlayer = [=](const double fraclastlayer_gauss, const int energy_bin) -> double
        {
            double fraclastlayer_norm = fraclastlayer_gauss;
            if (box_cox_correction_parameters.fraclast_norm_rms[energy_bin-1])
                fraclastlayer_norm /= box_cox_correction_parameters.fraclast_norm_rms[energy_bin-1];
            if (box_cox_correction_parameters.sumrms_norm_mean[energy_bin-1]>0)
                fraclastlayer_norm -= box_cox_correction_parameters.fraclast_norm_mean[energy_bin-1];
            else
                fraclastlayer_norm += abs(box_cox_correction_parameters.fraclast_norm_mean[energy_bin-1]);
            return fraclastlayer_norm;
        };

        auto normalize_xtrl = [=](const double xtrl_gauss, const int energy_bin) -> double
        {
            double xtrl_norm = xtrl_gauss;
            if (box_cox_correction_parameters.xtrl_norm_rms[energy_bin-1])
                xtrl_norm /= box_cox_correction_parameters.xtrl_norm_rms[energy_bin-1];
            if (box_cox_correction_parameters.xtrl_norm_mean[energy_bin-1]>0)
                xtrl_norm -= box_cox_correction_parameters.xtrl_norm_mean[energy_bin-1];
            else
                xtrl_norm += abs(box_cox_correction_parameters.xtrl_norm_mean[energy_bin-1]);
            return xtrl_norm;
        };

        my_vars.rms = normalize_rmslayer(my_vars.rms, get_energy_bin(corrected_energy_gev));
        my_vars.sumrms = normalize_sumrms(my_vars.sumrms, get_energy_bin(corrected_energy_gev));
        my_vars.fraclayer = normalize_fraclayer(my_vars.fraclayer, get_energy_bin(corrected_energy_gev));
        my_vars.fraclastlayer = normalize_fraclastlayer(my_vars.fraclastlayer, get_energy_bin(corrected_energy_gev));
        my_vars.xtrl = normalize_xtrl(my_vars.xtrl, get_energy_bin(corrected_energy_gev));

        return my_vars;
    }