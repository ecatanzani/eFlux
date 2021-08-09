#include "regularize.h"

#include "TF1.h"
#include "TVector3.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void load_tf1s(
    std::string file_path,
    std::vector<TF1> &sumrms_fitfunc,
    std::vector<TF1> &sumrms_fitfunc_err,
    std::vector<TF1> &flast_fitfunc,
    std::vector<TF1> &flast_fitfunc_err,
    const bool verbose)
{
    TFile *fitfile = TFile::Open(file_path.c_str(), "READ");
    if (!fitfile->IsOpen())
    {
        std::cerr << "\n\nError reading summary fit TTree [" << file_path << "]\n\n";
        exit(100);
    }
    else
        if (verbose)
            std::cout << "\nReading fitting summary [" << file_path << "]\n";

    // **** Fit function parameters
    const char* tree_name = "fit_summary";
    const char* func_form = "pol3";
    const int npars = 4;
    const double lvalue = 0;
    const double rvalue = 1; 
    // *****

    TTreeReader myReader(tree_name, fitfile);
    // Fit parameters
    TTreeReaderValue<std::vector<double>> flast_pars(myReader, "flast_pars");
    TTreeReaderValue<std::vector<double>> flast_err_pars(myReader, "flast_err_pars");
    TTreeReaderValue<std::vector<double>> sumrms_pars(myReader, "sumrms_pars");
    TTreeReaderValue<std::vector<double>> sumrms_err_pars(myReader, "sumrms_err_pars");
#if 0    
    // Fit ranges
    TTreeReaderValue<std::vector<double>> flast_fit_func_range(myReader, "flast_fit_func_range");
    TTreeReaderValue<std::vector<double>> flast_err_fit_func_range(myReader, "flast_err_fit_func_range");
    TTreeReaderValue<std::vector<double>> sumrms_fit_func_range(myReader, "sumrms_fit_func_range");
    TTreeReaderValue<std::vector<double>> sumrms_err_fit_func_range(myReader, "sumrms_err_fit_func_range");
#endif

    int bin_idx = 0;
    while (myReader.Next())
    {   
        sumrms_fitfunc[bin_idx] = TF1(
            (std::string(tree_name) + std::string("_sumrms_fitfunc_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);
        sumrms_fitfunc_err[bin_idx] = TF1(
            (std::string(tree_name) + std::string("_sumrms_fitfunc_err_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);
        flast_fitfunc[bin_idx] = TF1(
            (std::string(tree_name) + std::string("_flast_fitfunc_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);
        flast_fitfunc_err[bin_idx] = TF1(
            (std::string(tree_name) + std::string("_flast_fitfunc_err_") + std::to_string(bin_idx)).c_str(), 
            func_form, 
            lvalue, 
            rvalue);

        for (int par_idx=0; par_idx<npars; ++par_idx)
        {
            sumrms_fitfunc[bin_idx].SetParameter(par_idx, sumrms_pars->at(par_idx));
            sumrms_fitfunc_err[bin_idx].SetParameter(par_idx, sumrms_err_pars->at(par_idx));
            flast_fitfunc[bin_idx].SetParameter(par_idx, flast_pars->at(par_idx));
            flast_fitfunc_err[bin_idx].SetParameter(par_idx, flast_err_pars->at(par_idx));
        }
        ++bin_idx; 
    }

    fitfile->Close();
}