#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <iostream>

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

inline std::vector<std::shared_ptr<TH1D>> projection(TH2D* histo)
{
    std::vector<std::shared_ptr<TH1D>> projections (histo->GetNbinsX());
    for (int ibin=1; ibin<=histo->GetNbinsX(); ++ibin)
    {
        projections[ibin-1] = std::shared_ptr<TH1D>(histo->ProjectionY((std::string("projectionY_") + std::string(histo->GetName()) + std::string("_bin_") + std::to_string(ibin)).c_str(), ibin, ibin));
        TF1 fitfunc("fitfunc", "gaus", -0.2, 0.2);
        fitfunc.SetNpx(10000);
        projections[ibin-1]->SetLineWidth(3);
        projections[ibin-1]->Fit("fitfunc");
        fitfunc.SetLineWidth(3);
    }
    return projections;
}

void SliceNormalization(
    const char* full_histo_path, 
    const char* normalized_histo_path)
{
    auto slicenorm = [](TH2D* histo) { 
        for (int bidx=1; bidx<=histo->GetNbinsX(); ++bidx)
	{
		auto scale_factor = 1/histo->Integral(bidx, bidx, 1, histo->GetNbinsY());
		for (int bidy=1; bidy<=histo->GetNbinsY(); ++bidy)
		{
            		histo->SetBinContent(bidx, bidy, histo->GetBinContent(bidx, bidy)*scale_factor);
			histo->SetBinError(bidx, bidy, histo->GetBinError(bidx, bidy)*scale_factor);
        	}
	}
	return histo; };

    TFile *_input = TFile::Open(full_histo_path, "READ");
    if (_input->IsZombie())
    {
        std::cerr << "\n\nError opening input file: [" << full_histo_path << "]\n\n";
        exit(100);
    }

    // Simu histos
    auto h_energy_diff2D = slicenorm(static_cast<TH2D*>(_input->Get("Simu/h_energy_diff2D")));
    auto h_energy_diff2D_corr = slicenorm(static_cast<TH2D*>(_input->Get("Simu/h_energy_diff2D_corr")));
    auto h_energy_unfold = slicenorm(static_cast<TH2D*>(_input->Get("Simu/h_energy_unfold")));
    auto h_energy_unfold_corr = slicenorm(static_cast<TH2D*>(_input->Get("Simu/h_energy_unfold_corr")));

    auto h_energy_ps_diff2D = slicenorm(static_cast<TH2D*>(_input->Get("Preselection/Simu/h_energy_ps_diff2D")));
    auto h_energy_ps_diff2D_corr = slicenorm(static_cast<TH2D*>(_input->Get("Preselection/Simu/h_energy_ps_diff2D_corr")));
    auto h_energy_ps_unfold = slicenorm(static_cast<TH2D*>(_input->Get("Preselection/Simu/h_energy_ps_unfold")));
    auto h_energy_ps_unfold_corr = slicenorm(static_cast<TH2D*>(_input->Get("Preselection/Simu/h_energy_ps_unfold_corr")));

    auto projectionY_energy_diff2D = projection(h_energy_diff2D);
    auto projectionY_energy_diff2D_corr = projection(h_energy_diff2D_corr);

    auto projectionY_ps_energy_diff2D = projection(h_energy_ps_diff2D);
    auto projectionY_ps_energy_diff2D_corr = projection(h_energy_ps_diff2D_corr);   

    h_energy_diff2D->SetDirectory(0);
    h_energy_diff2D_corr->SetDirectory(0);
    h_energy_unfold->SetDirectory(0);
    h_energy_unfold_corr->SetDirectory(0);

    h_energy_ps_diff2D->SetDirectory(0);
    h_energy_ps_diff2D_corr->SetDirectory(0);
    h_energy_ps_unfold->SetDirectory(0);
    h_energy_ps_unfold_corr->SetDirectory(0);

    // XTRL histo
    auto h_xtrl = slicenorm(static_cast<TH2D*>(_input->Get("BGO/h_xtrl")));
    auto h_xtrl_ps = slicenorm(static_cast<TH2D*>(_input->Get("Preselection/BGO/h_xtrl_ps")));

    h_xtrl->SetDirectory(0);
    h_xtrl_ps->SetDirectory(0);
    
    TFile *_output = TFile::Open(normalized_histo_path, "RECREATE");
    if (_output->IsZombie())
    {
        std::cerr << "\n\nError writing output file: [" << normalized_histo_path << "]\n\n";
        exit(100);
    }

    _output->mkdir("Simu");
    _output->cd("Simu");

    h_energy_diff2D->Write();
    h_energy_diff2D_corr->Write();
    h_energy_unfold->Write();
    h_energy_unfold_corr->Write();

    _output->mkdir("Simu/Projections");
    _output->cd("Simu/Projections");
    for (auto& _elm : projectionY_energy_diff2D)
        _elm->Write();
    for (auto& _elm : projectionY_energy_diff2D_corr)
        _elm->Write();

    _output->mkdir("xtrl");
    _output->cd("xtrl");

    h_xtrl->Write();

    _output->mkdir("Preselection/Simu");
    _output->cd("Preselection/Simu");

    h_energy_ps_diff2D->Write();
    h_energy_ps_diff2D_corr->Write();
    h_energy_ps_unfold->Write();
    h_energy_ps_unfold_corr->Write();

    _output->mkdir("Preselection/Simu/Projections");
    _output->cd("Preselection/Simu/Projections");

    for (auto& _elm : projectionY_ps_energy_diff2D)
        _elm->Write();
    for (auto& _elm : projectionY_ps_energy_diff2D_corr)
        _elm->Write();    

    _output->mkdir("Preselection/xtrl");
    _output->cd("Preselection/xtrl");

    h_xtrl_ps->Write();

    _output->Close();

}

inline double get_left_index(TProfile* profile)
{
    double firstx = 0;
    for (int idx=1; idx<=profile->GetNbinsX(); ++idx)
        if (profile->GetBinContent(idx))
        {
            firstx = profile->GetBinCenter(idx);
            break;
        }
    return firstx;
}

void FitSummary(
    std::vector<std::shared_ptr<TF1>> fit_func, 
    const char* output_file_name,
    const char* output_tree_name)
{
    TFile* fit_summary = TFile::Open(output_file_name, "RECREATE");
    if (!fit_summary->IsOpen())
    {
        std::cerr << "\n\nError writing summary fit TFile [" << output_file_name << "]\n\n";
        exit(100);
    }

    TTree fit_tree(output_tree_name, "Fit Summary Tree");

    std::vector<double> pars (fit_func[0]->GetNpar());
    double chi2, chi2_ndf;
    int ndof;

    fit_tree.Branch("pars", &pars);
    fit_tree.Branch("chi2", &chi2, "chi2/D");
    fit_tree.Branch("chi2_ndf", &chi2_ndf, "chi2_ndf/D");
    fit_tree.Branch("ndof", &ndof, "ndof/I");

    for (unsigned int idx=0; idx<fit_func.size(); ++idx)
    {
        for (int idx_par = 0; idx_par < fit_func[idx]->GetNpar(); ++idx_par)
            pars[idx_par] =  fit_func[idx]->GetParameter(idx_par);
            
        chi2 = fit_func[idx]->GetChisquare();
        ndof = fit_func[idx]->GetNDF();
        chi2_ndf = chi2/ndof;

        fit_tree.Fill();
    }

    fit_tree.Write();

    fit_summary->Close();
}

void FullFitSummary(
    std::vector<std::shared_ptr<TF1>> flast_profile_fit_func,
    std::vector<std::shared_ptr<TF1>> flast_err_fit_func,
    std::vector<std::shared_ptr<TF1>> sumrms_profile_fit_func,
    std::vector<std::shared_ptr<TF1>> sumrms_err_fit_func,
    const char* output_file_name,
    const char* output_tree_name)
{
    TFile* fit_summary = TFile::Open(output_file_name, "RECREATE");
    if (!fit_summary->IsOpen())
    {
        std::cerr << "\n\nError writing summary fit TFile [" << output_file_name << "]\n\n";
        exit(100);
    }

    TTree fit_tree(output_tree_name, "Fit Summary Tree");

    // Get number of fit parameters
    const auto npars = flast_profile_fit_func[0]->GetNpar();
    // Get number of energy bins
    const auto nbins = flast_profile_fit_func.size();
    // vector to store fit parameters
    std::vector<double> flast_pars (npars, 0);
    std::vector<double> flast_err_pars (npars, 0);
    std::vector<double> sumrms_pars (npars, 0);
    std::vector<double> sumrms_err_pars (npars, 0);
    // Chi2
    double flast_chi2;
    double flast_chi2_ndf;
    double flast_err_chi2;
    double flast_err_chi2_ndf;
    double sumrms_chi2;
    double sumrms_chi2_ndf;
    double sumrms_err_chi2;
    double sumrms_err_chi2_ndf;
    // dof
    int flast_ndof;
    int flast_err_ndof;
    int sumrms_ndof;
    int sumrms_err_ndof;
    // Functions range
    std::vector<double> sumrms_fit_func_range(2, 0);
    std::vector<double> sumrms_err_fit_func_range(2, 0);
    std::vector<double> flast_fit_func_range(2, 0);
    std::vector<double> flast_err_fit_func_range(2, 0);

    // Fit parameters
    fit_tree.Branch("flast_pars", &flast_pars);
    fit_tree.Branch("flast_err_pars", &flast_err_pars);
    fit_tree.Branch("sumrms_pars", &sumrms_pars);
    fit_tree.Branch("sumrms_err_pars", &sumrms_err_pars);

    // Chi2
    fit_tree.Branch("flast_chi2", &flast_chi2, "flast_chi2/D");
    fit_tree.Branch("flast_chi2_ndf", &flast_chi2_ndf, "flast_chi2_ndf/D");
    fit_tree.Branch("flast_err_chi2", &flast_err_chi2, "flast_err_chi2/D");
    fit_tree.Branch("flast_err_chi2_ndf", &flast_err_chi2_ndf, "flast_err_chi2_ndf/D");

    fit_tree.Branch("sumrms_chi2", &sumrms_chi2, "sumrms_chi2/D");
    fit_tree.Branch("sumrms_chi2_ndf", &sumrms_chi2_ndf, "sumrms_chi2_ndf/D");
    fit_tree.Branch("sumrms_err_chi2", &sumrms_err_chi2, "sumrms_err_chi2/D");
    fit_tree.Branch("sumrms_err_chi2_ndf", &sumrms_err_chi2_ndf, "sumrms_err_chi2_ndf/D");

    // dof    
    fit_tree.Branch("flast_ndof", &flast_ndof, "flast_ndof/I");
    fit_tree.Branch("flast_err_ndof", &flast_err_ndof, "flast_err_ndof/I");
    fit_tree.Branch("sumrms_ndof", &sumrms_ndof, "sumrms_ndof/I");
    fit_tree.Branch("sumrms_err_ndof", &sumrms_err_ndof, "sumrms_err_ndof/I");

    // Functions range    
    fit_tree.Branch("flast_fit_func_range", &flast_fit_func_range);
    fit_tree.Branch("flast_err_fit_func_range", &flast_err_fit_func_range);
    fit_tree.Branch("sumrms_fit_func_range", &sumrms_fit_func_range);
    fit_tree.Branch("sumrms_err_fit_func_range", &sumrms_err_fit_func_range);


    for (unsigned int idx=0; idx<nbins; ++idx)
    {
        for (int idx_par = 0; idx_par < npars; ++idx_par)
        {
            flast_pars[idx_par] =  flast_profile_fit_func[idx]->GetParameter(idx_par);
            flast_err_pars[idx_par] =  flast_err_fit_func[idx]->GetParameter(idx_par);
            sumrms_pars[idx_par] =  sumrms_profile_fit_func[idx]->GetParameter(idx_par);
            sumrms_err_pars[idx_par] =  sumrms_err_fit_func[idx]->GetParameter(idx_par);
        }

        flast_fit_func_range[0] = flast_profile_fit_func[idx]->GetXmin();
        flast_fit_func_range[1] = flast_profile_fit_func[idx]->GetXmax();
        flast_err_fit_func_range[0] = flast_err_fit_func[idx]->GetXmin();
        flast_err_fit_func_range[1] = flast_err_fit_func[idx]->GetXmax();
        sumrms_fit_func_range[0] = sumrms_profile_fit_func[idx]->GetXmin();
        sumrms_fit_func_range[1] = sumrms_profile_fit_func[idx]->GetXmax();
        sumrms_err_fit_func_range[0] = sumrms_err_fit_func[idx]->GetXmin();
        sumrms_err_fit_func_range[1] = sumrms_err_fit_func[idx]->GetXmax();


        flast_chi2 = flast_profile_fit_func[idx]->GetChisquare();
        flast_err_chi2 = flast_err_fit_func[idx]->GetChisquare();
        sumrms_chi2 = sumrms_profile_fit_func[idx]->GetChisquare();
        sumrms_err_chi2 = sumrms_err_fit_func[idx]->GetChisquare();

        flast_ndof = flast_profile_fit_func[idx]->GetNDF();
        flast_err_ndof = flast_err_fit_func[idx]->GetNDF();
        sumrms_ndof = sumrms_profile_fit_func[idx]->GetNDF();
        sumrms_err_ndof = sumrms_err_fit_func[idx]->GetNDF();

        flast_chi2_ndf = flast_chi2/flast_ndof;
        flast_err_chi2_ndf = flast_err_chi2/flast_err_ndof;
        sumrms_chi2_ndf = sumrms_chi2/sumrms_ndof;
        sumrms_err_chi2_ndf = sumrms_err_chi2/sumrms_err_ndof;

        fit_tree.Fill();
    }

    fit_tree.Write();

    fit_summary->Close();
}


void FlastCosineProfile(
    const char* full_histo_path, 
    const char* output,
    const int nbins_energy=50,
    const char* profile_err_opt = "s",
    const bool full_profiles_canvas = true,
    const bool verbose = true)
{
    std::vector<TH2D*> flast_cosine (nbins_energy);
    std::vector<TProfile*> flast_cosine_profile (nbins_energy);
    std::vector<std::shared_ptr<TH1D>> flast_cosine_err (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_fit_func_p2 (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_fit_func_p3 (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_fit_func_p4 (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_err_fit_func_p2 (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_err_fit_func_p3 (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_err_fit_func_p4 (nbins_energy);

    if (verbose)
        std::cout << "\nReading input ROOT file [" << full_histo_path << "]\n";
    TFile *input_file = TFile::Open(full_histo_path, "READ");
    if (input_file->IsZombie())
    {
        std::cerr << "\n\nError opening input file: [" << full_histo_path << "]\n\n";
        exit(100);
    }
    
    if (verbose)
        std::cout << "\nFitting... \n";

    for (int idx=0; idx<nbins_energy; ++idx)
    {   
        // Read 2D histo
        std::string hname = "Preselection/BGO/energybin_" + to_string(idx+1) + "/h_BGOrec_ps_ratio_last_cosine2D_fdr_bin_" + std::to_string(idx+1);

        flast_cosine[idx] = static_cast<TH2D*>(input_file->Get(hname.c_str()));
        flast_cosine[idx]->SetDirectory(0);
        
        // Profile
        flast_cosine_profile[idx] = static_cast<TProfile*>(flast_cosine[idx]->ProfileX(
            (std::string(flast_cosine[idx]->GetName()) + "_profileX").c_str() ,0 , flast_cosine[idx]->GetNbinsY(), profile_err_opt));
        flast_cosine_profile[idx]->SetDirectory(0);

        flast_cosine_err[idx] = std::make_shared<TH1D>(
            (std::string("profileX_err_") + std::to_string(idx+1)).c_str(), 
            "profileX Error; cos(#theta); RMS", 
            flast_cosine_profile[idx]->GetXaxis()->GetNbins(), 0, 1);
       
        for (int b_idx=1; b_idx<=flast_cosine_err[idx]->GetNbinsX(); ++b_idx)
        {
            flast_cosine_err[idx]->SetBinContent(b_idx, flast_cosine_profile[idx]->GetBinError(b_idx));
            flast_cosine_err[idx]->SetBinError(b_idx, sqrt(flast_cosine_profile[idx]->GetBinError(b_idx)));
        }

        // TF1s
        auto lidx = get_left_index(flast_cosine_profile[idx]);
        profile_fit_func_p2[idx] = std::make_shared<TF1>((std::string("profile_fitfunc_p2_") + std::to_string(idx+1)).c_str(), "pol2", lidx, 1);
        profile_fit_func_p3[idx] = std::make_shared<TF1>((std::string("profile_fitfunc_p3_") + std::to_string(idx+1)).c_str(), "pol3", lidx, 1);
        profile_fit_func_p4[idx] = std::make_shared<TF1>((std::string("profile_fitfunc_p4_") + std::to_string(idx+1)).c_str(), "pol4", lidx, 1); 

        profile_err_fit_func_p2[idx] = std::make_shared<TF1>((std::string("profile_err_fitfunc_p2_") + std::to_string(idx+1)).c_str(), "pol2", lidx, 1);
        profile_err_fit_func_p3[idx] = std::make_shared<TF1>((std::string("profile_err_fitfunc_p3_") + std::to_string(idx+1)).c_str(), "pol3", lidx, 1);
        profile_err_fit_func_p4[idx] = std::make_shared<TF1>((std::string("profile_err_fitfunc_p4_") + std::to_string(idx+1)).c_str(), "pol4", lidx, 1); 
        
        profile_fit_func_p2[idx]->SetNpx(1000);
        profile_fit_func_p3[idx]->SetNpx(1000);
        profile_fit_func_p4[idx]->SetNpx(1000);
        profile_err_fit_func_p2[idx]->SetNpx(1000);
        profile_err_fit_func_p3[idx]->SetNpx(1000);
        profile_err_fit_func_p4[idx]->SetNpx(1000);

        profile_fit_func_p2[idx]->SetLineColor(kGreen);
        profile_fit_func_p3[idx]->SetLineColor(kCyan);
        profile_fit_func_p4[idx]->SetLineColor(kRed);
        profile_err_fit_func_p2[idx]->SetLineColor(kGreen);
        profile_err_fit_func_p3[idx]->SetLineColor(kCyan);
        profile_err_fit_func_p4[idx]->SetLineColor(kRed);
        
        // Fitting
        flast_cosine_profile[idx]->Fit(profile_fit_func_p2[idx]->GetName(), "WQRN");
        flast_cosine_profile[idx]->Fit(profile_fit_func_p3[idx]->GetName(), "WQRN");
        flast_cosine_profile[idx]->Fit(profile_fit_func_p4[idx]->GetName(), "WQRN");

        flast_cosine_err[idx]->Fit(profile_err_fit_func_p2[idx]->GetName(), "WQRN");
        flast_cosine_err[idx]->Fit(profile_err_fit_func_p3[idx]->GetName(), "WQRN");
        flast_cosine_err[idx]->Fit(profile_err_fit_func_p4[idx]->GetName(), "WQRN");
    }
    input_file->Close();
    
    if (verbose)
        std::cout << "\nWriting output ROOT file [" << output << "]\n";
    TFile *output_file = TFile::Open(output, "RECREATE");
    if (output_file->IsZombie())
    {
        std::cerr << "\n\nError writing output file: [" << output << "]\n\n";
        exit(100);
    }

    for (int idx=0; idx<nbins_energy; ++idx)
    {    
        output_file->mkdir((std::string("energybin_") + std::to_string(idx+1)).c_str());
        output_file->cd((std::string("energybin_") + std::to_string(idx+1)).c_str());
        flast_cosine[idx]->Write();
        flast_cosine_profile[idx]->Write();
        flast_cosine_err[idx]->Write();
        profile_fit_func_p2[idx]->Write();
        profile_fit_func_p3[idx]->Write();
        profile_fit_func_p4[idx]->Write();
        profile_err_fit_func_p2[idx]->Write();
        profile_err_fit_func_p3[idx]->Write();
        profile_err_fit_func_p4[idx]->Write();
    }
    
    if (full_profiles_canvas)
    {
        output_file->cd();
        std::unique_ptr<TCanvas> c_full_profiles = std::make_unique<TCanvas>("c_full_profiles", "Full Profiles", 1500, 950);
        std::unique_ptr<TCanvas> c_full_profiles_err = std::make_unique<TCanvas>("c_full_profiles_err", "Full Profiles Errors", 1500, 950);
        c_full_profiles->Divide(5, 10);
        c_full_profiles_err->Divide(5, 10);
        
        for (int idx=0; idx<nbins_energy; ++idx)
        {
            // Set ranges
            flast_cosine_profile[idx]->GetXaxis()->SetRangeUser(0.6, 1);
            flast_cosine_err[idx]->GetXaxis()->SetRangeUser(0.6, 1);
            
            // Profile canvas
            c_full_profiles->cd(idx+1);
            gStyle->SetOptStat(0);
            flast_cosine_profile[idx]->Draw();
            profile_fit_func_p2[idx]->Draw("same");
            profile_fit_func_p3[idx]->Draw("same");
            profile_fit_func_p4[idx]->Draw("same");

            // Profile err canvas
            c_full_profiles_err->cd(idx+1);
            gStyle->SetOptStat(0);
            flast_cosine_err[idx]->Draw();
            profile_err_fit_func_p2[idx]->Draw("same");
            profile_err_fit_func_p3[idx]->Draw("same");
            profile_err_fit_func_p4[idx]->Draw("same");
        }

        c_full_profiles->Write();
        c_full_profiles_err->Write();
    }

    output_file->Close();
    
    // Write fit summary
    if (verbose)
        std::cout << "\nWriting pol2 fit summary [fit_pol2_summary.root]\n";
    FitSummary(profile_fit_func_p2, "fit_pol2_summary.root", "fit_pol2_tree");
    if (verbose)
        std::cout << "\nWriting pol3 fit summary [fit_pol3_summary.root]\n";
    FitSummary(profile_fit_func_p3, "fit_pol3_summary.root", "fit_pol3_tree");
    if (verbose)
        std::cout << "\nWriting pol2 fit summary [fit_pol4_summary.root]\n";
    FitSummary(profile_fit_func_p4, "fit_pol4_summary.root", "fit_pol4_tree");
    if (verbose)
        std::cout << "\nWriting pol2 error fit summary [fit_pol2_err_summary.root]\n";
    FitSummary(profile_err_fit_func_p2, "fit_pol2_err_summary.root", "fit_pol2_err_tree");
    if (verbose)
        std::cout << "\nWriting pol3 error fit summary [fit_pol3_err_summary.root]\n";
    FitSummary(profile_err_fit_func_p3, "fit_pol3_err_summary.root", "fit_pol3_err_tree");
    if (verbose)
        std::cout << "\nWriting pol2 error fit summary [fit_pol4_err_summary.root]\n";
    FitSummary(profile_err_fit_func_p4, "fit_pol4_err_summary.root", "fit_pol4_err_tree");
    
}


void sumRMSCosineProfile(
     const char* full_histo_path, 
    const char* output,
    const int nbins_energy=50,
    const char* profile_err_opt = "s",
    const bool full_profiles_canvas = true,
    const bool verbose = true)
{
    std::vector<TH2D*> sumrms_cosine (nbins_energy);
    std::vector<TProfile*> sumrms_cosine_profile (nbins_energy);
    std::vector<std::shared_ptr<TH1D>> sumrms_cosine_err (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_fit_func_p2 (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_fit_func_p3 (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_fit_func_p4 (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_err_fit_func_p2 (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_err_fit_func_p3 (nbins_energy);
    std::vector<std::shared_ptr<TF1>> profile_err_fit_func_p4 (nbins_energy);

    if (verbose)
        std::cout << "\nReading input ROOT file [" << full_histo_path << "]\n";
    TFile *input_file = TFile::Open(full_histo_path, "READ");
    if (input_file->IsZombie())
    {
        std::cerr << "\n\nError opening input file: [" << full_histo_path << "]\n\n";
        exit(100);
    }
    
    if (verbose)
        std::cout << "\nFitting... \n";

    for (int idx=0; idx<nbins_energy; ++idx)
    {   
        // Read 2D histo
        std::string hname = "Preselection/BGO/energybin_" + to_string(idx+1) + "/h_BGOrec_ps_sumRms_cosine2D_bin_" + std::to_string(idx+1);
        std::string hname_norm = "h_BGOrec_ps_sumRms_cosine2D_norm_bin_" + std::to_string(idx+1);
        sumrms_cosine[idx] = static_cast<TH2D*>(input_file->Get(hname.c_str()));
        sumrms_cosine[idx]->SetDirectory(0);
        
        // Profile
        sumrms_cosine_profile[idx] = static_cast<TProfile*>(sumrms_cosine[idx]->ProfileX(
            (std::string(sumrms_cosine[idx]->GetName()) + "_profileX").c_str() ,0 , sumrms_cosine[idx]->GetYaxis()->GetNbins()));
        sumrms_cosine_profile[idx]->SetDirectory(0);

        sumrms_cosine_err[idx] = std::make_shared<TH1D>(
            (std::string("profileX_err_") + std::to_string(idx+1)).c_str(), 
            "profileX Error; cos(#theta); RMS", 
            sumrms_cosine_profile[idx]->GetXaxis()->GetNbins(), 0, 1);

        for (int b_idx=1; b_idx<=sumrms_cosine_err[idx]->GetNbinsX(); ++b_idx)
        {
            sumrms_cosine_err[idx]->SetBinContent(b_idx, sumrms_cosine_profile[idx]->GetBinError(b_idx));
            sumrms_cosine_err[idx]->SetBinError(b_idx, sqrt(sumrms_cosine_profile[idx]->GetBinError(b_idx)));
        }

        // TF1s
        auto lidx = get_left_index(sumrms_cosine_profile[idx]);
        profile_fit_func_p2[idx] = std::make_shared<TF1>((std::string("profile_fitfunc_p2_") + std::to_string(idx+1)).c_str(), "pol2", lidx, 1);
        profile_fit_func_p3[idx] = std::make_shared<TF1>((std::string("profile_fitfunc_p3_") + std::to_string(idx+1)).c_str(), "pol3", lidx, 1);
        profile_fit_func_p4[idx] = std::make_shared<TF1>((std::string("profile_fitfunc_p4_") + std::to_string(idx+1)).c_str(), "pol4", lidx, 1); 

        profile_err_fit_func_p2[idx] = std::make_shared<TF1>((std::string("profile_err_fitfunc_p2_") + std::to_string(idx+1)).c_str(), "pol2", lidx, 1);
        profile_err_fit_func_p3[idx] = std::make_shared<TF1>((std::string("profile_err_fitfunc_p3_") + std::to_string(idx+1)).c_str(), "pol3", lidx, 1);
        profile_err_fit_func_p4[idx] = std::make_shared<TF1>((std::string("profile_err_fitfunc_p4_") + std::to_string(idx+1)).c_str(), "pol4", lidx, 1);

        profile_fit_func_p2[idx]->SetNpx(1000);
        profile_fit_func_p3[idx]->SetNpx(1000);
        profile_fit_func_p4[idx]->SetNpx(1000);
        profile_err_fit_func_p2[idx]->SetNpx(1000);
        profile_err_fit_func_p3[idx]->SetNpx(1000);
        profile_err_fit_func_p4[idx]->SetNpx(1000);

        profile_fit_func_p2[idx]->SetLineColor(kGreen);
        profile_fit_func_p3[idx]->SetLineColor(kCyan);
        profile_fit_func_p4[idx]->SetLineColor(kRed);
        profile_err_fit_func_p2[idx]->SetLineColor(kGreen);
        profile_err_fit_func_p3[idx]->SetLineColor(kCyan);
        profile_err_fit_func_p4[idx]->SetLineColor(kRed);

        // Fitting
        sumrms_cosine_profile[idx]->Fit(profile_fit_func_p2[idx]->GetName(), "WQRN");
        sumrms_cosine_profile[idx]->Fit(profile_fit_func_p3[idx]->GetName(), "WQRN");
        sumrms_cosine_profile[idx]->Fit(profile_fit_func_p4[idx]->GetName(), "WQRN");

        sumrms_cosine_err[idx]->Fit(profile_err_fit_func_p2[idx]->GetName(), "WQRN");
        sumrms_cosine_err[idx]->Fit(profile_err_fit_func_p3[idx]->GetName(), "WQRN");
        sumrms_cosine_err[idx]->Fit(profile_err_fit_func_p4[idx]->GetName(), "WQRN");
    }
    input_file->Close();

    if (verbose)
        std::cout << "\nWriting output ROOT file [" << output << "]\n";
    TFile *output_file = TFile::Open(output, "RECREATE");
    if (output_file->IsZombie())
    {
        std::cerr << "\n\nError writing output file: [" << output << "]\n\n";
        exit(100);
    }

    for (int idx=0; idx<nbins_energy; ++idx)
    {    
        output_file->mkdir((std::string("energybin_") + std::to_string(idx+1)).c_str());
        output_file->cd((std::string("energybin_") + std::to_string(idx+1)).c_str());
        sumrms_cosine[idx]->Write();
        sumrms_cosine_profile[idx]->Write();
        sumrms_cosine_err[idx]->Write();
        profile_fit_func_p2[idx]->Write();
        profile_fit_func_p3[idx]->Write();
        profile_fit_func_p4[idx]->Write();
        profile_err_fit_func_p2[idx]->Write();
        profile_err_fit_func_p3[idx]->Write();
        profile_err_fit_func_p4[idx]->Write();
    }

    if (full_profiles_canvas)
    {
        output_file->cd();
        std::unique_ptr<TCanvas> c_full_profiles = std::make_unique<TCanvas>("c_full_profiles", "Full Profiles", 1500, 950);
        std::unique_ptr<TCanvas> c_full_profiles_err = std::make_unique<TCanvas>("c_full_profiles_err", "Full Profiles Errors", 1500, 950);
        c_full_profiles->Divide(5, 10);
        c_full_profiles_err->Divide(5, 10);

        for (int idx=0; idx<nbins_energy; ++idx)
        {
            // Set ranges
            sumrms_cosine_profile[idx]->GetXaxis()->SetRangeUser(0.6, 1);
            sumrms_cosine_err[idx]->GetXaxis()->SetRangeUser(0.6, 1);
            
            // Profile canvas
            c_full_profiles->cd(idx+1);
            gStyle->SetOptStat(0);
            sumrms_cosine_profile[idx]->Draw();
            profile_fit_func_p2[idx]->Draw("same");
            profile_fit_func_p3[idx]->Draw("same");
            profile_fit_func_p4[idx]->Draw("same");

            // Profile err canvas
            c_full_profiles_err->cd(idx+1);
            gStyle->SetOptStat(0);
            sumrms_cosine_err[idx]->Draw();
            profile_err_fit_func_p2[idx]->Draw("same");
            profile_err_fit_func_p3[idx]->Draw("same");
            profile_err_fit_func_p4[idx]->Draw("same");
        }

        c_full_profiles->Write();
        c_full_profiles_err->Write();
    }

    output_file->Close();

    // Write fit summary
    if (verbose)
        std::cout << "\nWriting pol2 fit summary [fit_pol2_summary.root]\n";
    FitSummary(profile_fit_func_p2, "fit_pol2_summary.root", "fit_pol2_tree");
    if (verbose)
        std::cout << "\nWriting pol3 fit summary [fit_pol3_summary.root]\n";
    FitSummary(profile_fit_func_p3, "fit_pol3_summary.root", "fit_pol3_tree");
    if (verbose)
        std::cout << "\nWriting pol2 fit summary [fit_pol4_summary.root]\n";
    FitSummary(profile_fit_func_p4, "fit_pol4_summary.root", "fit_pol4_tree");
    if (verbose)
        std::cout << "\nWriting pol2 error fit summary [fit_pol2_err_summary.root]\n";
    FitSummary(profile_err_fit_func_p2, "fit_pol2_err_summary.root", "fit_pol2_err_tree");
    if (verbose)
        std::cout << "\nWriting pol3 error fit summary [fit_pol3_err_summary.root]\n";
    FitSummary(profile_err_fit_func_p3, "fit_pol3_err_summary.root", "fit_pol3_err_tree");
    if (verbose)
        std::cout << "\nWriting pol2 error fit summary [fit_pol4_err_summary.root]\n";
    FitSummary(profile_err_fit_func_p4, "fit_pol4_err_summary.root", "fit_pol4_err_tree");
}

void BuildVariablesProfiles(
    const char* full_histo_path, 
    const char* output,
    const char* fit_summary_output,
    const int nbins_energy=50,
    const char* profile_err_opt = "s",
    const bool full_profiles_canvas = true,
    const bool verbose = true)
{
    std::vector<TH2D*> sumrms_cosine (nbins_energy);
    std::vector<TProfile*> sumrms_cosine_profile (nbins_energy);
    std::vector<std::shared_ptr<TH1D>> sumrms_cosine_err (nbins_energy);
    std::vector<std::shared_ptr<TF1>> sumrms_profile_fit_func (nbins_energy);
    std::vector<std::shared_ptr<TF1>> sumrms_profile_err_fit_func (nbins_energy);

    std::vector<TH2D*> flast_cosine (nbins_energy);
    std::vector<TProfile*> flast_cosine_profile (nbins_energy);
    std::vector<std::shared_ptr<TH1D>> flast_cosine_err (nbins_energy);
    std::vector<std::shared_ptr<TF1>> flast_profile_fit_func (nbins_energy);
    std::vector<std::shared_ptr<TF1>> flast_profile_err_fit_func (nbins_energy);

    if (verbose)
        std::cout << "\nReading input ROOT file [" << full_histo_path << "]\n";
    TFile *input_file = TFile::Open(full_histo_path, "READ");
    if (input_file->IsZombie())
    {
        std::cerr << "\n\nError opening input file: [" << full_histo_path << "]\n\n";
        exit(100);
    }

    if (verbose)
        std::cout << "\nFitting... \n";

    for (int idx=0; idx<nbins_energy; ++idx)
    {   
        // Read 2D histo
        std::string hname = "Preselection/BGO/energybin_" + to_string(idx+1) + "/h_BGOrec_ps_ratio_last_cosine2D_fdr_bin_" + std::to_string(idx+1);
        flast_cosine[idx] = static_cast<TH2D*>(input_file->Get(hname.c_str()));
        hname = "Preselection/BGO/energybin_" + to_string(idx+1) + "/h_BGOrec_ps_sumRms_cosine2D_bin_" + std::to_string(idx+1);
        sumrms_cosine[idx] = static_cast<TH2D*>(input_file->Get(hname.c_str()));
        flast_cosine[idx]->SetDirectory(0);
        sumrms_cosine[idx]->SetDirectory(0);

        // Build profiles
        flast_cosine_profile[idx] = static_cast<TProfile*>(flast_cosine[idx]->ProfileX(
            (std::string(flast_cosine[idx]->GetName()) + "_profileX").c_str() ,0 , flast_cosine[idx]->GetNbinsY(), profile_err_opt));
        sumrms_cosine_profile[idx] = static_cast<TProfile*>(sumrms_cosine[idx]->ProfileX(
            (std::string(sumrms_cosine[idx]->GetName()) + "_profileX").c_str() ,0 , sumrms_cosine[idx]->GetNbinsY(), profile_err_opt));
        flast_cosine_profile[idx]->SetDirectory(0);
        sumrms_cosine_profile[idx]->SetDirectory(0);

        flast_cosine_err[idx] = std::make_shared<TH1D>(
            (std::string(flast_cosine[idx]->GetName()) + std::string("_profileX_err_") + std::to_string(idx+1)).c_str(), 
            "profileX Error; cos(#theta); RMS", 
            flast_cosine_profile[idx]->GetNbinsX(), 0, 1);
        sumrms_cosine_err[idx] = std::make_shared<TH1D>(
            (std::string(sumrms_cosine[idx]->GetName()) + std::string("_profileX_err_") + std::to_string(idx+1)).c_str(), 
            "profileX Error; cos(#theta); RMS", 
            sumrms_cosine_profile[idx]->GetNbinsX(), 0, 1);
       
        for (int b_idx=1; b_idx<=flast_cosine_err[idx]->GetNbinsX(); ++b_idx)
        {
            flast_cosine_err[idx]->SetBinContent(b_idx, flast_cosine_profile[idx]->GetBinError(b_idx));
            auto tmp_profileY = flast_cosine[idx]->ProfileY((std::string(flast_cosine[idx]->GetName()) + "_profileY").c_str() ,b_idx , b_idx, profile_err_opt);
            flast_cosine_err[idx]->SetBinError(b_idx, tmp_profileY->GetRMSError());
            
            sumrms_cosine_err[idx]->SetBinContent(b_idx, sumrms_cosine_profile[idx]->GetBinError(b_idx));
            tmp_profileY = sumrms_cosine[idx]->ProfileY((std::string(sumrms_cosine[idx]->GetName()) + "_profileY").c_str() ,b_idx , b_idx, profile_err_opt);
            sumrms_cosine_err[idx]->SetBinError(b_idx, tmp_profileY->GetRMSError());
        }

        // TF1s
        auto flast_lidx = get_left_index(flast_cosine_profile[idx]);
        auto sumrms_lidx = get_left_index(flast_cosine_profile[idx]);
        flast_profile_fit_func[idx] = std::make_shared<TF1>((std::string("flast_profile_fitfunc_") + std::to_string(idx+1)).c_str(), "pol3", flast_lidx, 1);
        flast_profile_err_fit_func[idx] = std::make_shared<TF1>((std::string("flast_profile_err_fitfunc_") + std::to_string(idx+1)).c_str(), "pol3", flast_lidx, 1);
        sumrms_profile_fit_func[idx] = std::make_shared<TF1>((std::string("sumrms_profile_fitfunc_") + std::to_string(idx+1)).c_str(), "pol3", sumrms_lidx, 1);
        sumrms_profile_err_fit_func[idx] = std::make_shared<TF1>((std::string("sumrms_profile_err_fitfunc_") + std::to_string(idx+1)).c_str(), "pol3", sumrms_lidx, 1);

        flast_profile_fit_func[idx]->SetNpx(1000);
        flast_profile_err_fit_func[idx]->SetNpx(1000);
        sumrms_profile_fit_func[idx]->SetNpx(1000);
        sumrms_profile_err_fit_func[idx]->SetNpx(1000);
        
        // Fitting
        flast_cosine_profile[idx]->Fit(flast_profile_fit_func[idx]->GetName(), "WQRN");
        flast_cosine_err[idx]->Fit(flast_profile_err_fit_func[idx]->GetName(), "WQRN");
        sumrms_cosine_profile[idx]->Fit(sumrms_profile_fit_func[idx]->GetName(), "WQRN");
        sumrms_cosine_err[idx]->Fit(sumrms_profile_err_fit_func[idx]->GetName(), "WQRN");
    }
    
    input_file->Close();

    if (verbose)
        std::cout << "\nWriting output ROOT file [" << output << "]\n";
    TFile *output_file = TFile::Open(output, "RECREATE");
    if (output_file->IsZombie())
    {
        std::cerr << "\n\nError writing output file: [" << output << "]\n\n";
        exit(100);
    }    

    for (int idx=0; idx<nbins_energy; ++idx)
    {    
        output_file->mkdir((std::string("energybin_") + std::to_string(idx+1)).c_str());
        output_file->cd((std::string("energybin_") + std::to_string(idx+1)).c_str());
        // Save flast results
        flast_cosine[idx]->Write();
        flast_cosine_profile[idx]->Write();
        flast_cosine_err[idx]->Write();
        flast_profile_fit_func[idx]->Write();
        flast_profile_err_fit_func[idx]->Write();
        // Save sumrms results
        sumrms_cosine[idx]->Write();
        sumrms_cosine_profile[idx]->Write();
        sumrms_cosine_err[idx]->Write();
        sumrms_profile_fit_func[idx]->Write();
        sumrms_profile_err_fit_func[idx]->Write();
    }

    if (full_profiles_canvas)
    {
        output_file->cd();
        double lidx = 0.6;
        std::unique_ptr<TCanvas> flast_full_profiles = std::make_unique<TCanvas>("flast_full_profiles", "Profiles", 1500, 950);
        std::unique_ptr<TCanvas> flast_full_profiles_err = std::make_unique<TCanvas>("flast_full_profiles_err", "Profiles Errors", 1500, 950);
        std::unique_ptr<TCanvas> sumrms_full_profiles = std::make_unique<TCanvas>("sumrms_full_profiles", "Profiles", 1500, 950);
        std::unique_ptr<TCanvas> sumrms_full_profiles_err = std::make_unique<TCanvas>("sumrms_full_profiles_err", "Profiles Errors", 1500, 950);
        flast_full_profiles->Divide(5, 10);
        flast_full_profiles_err->Divide(5, 10);
        sumrms_full_profiles->Divide(5, 10);
        sumrms_full_profiles_err->Divide(5, 10);

        for (int idx=0; idx<nbins_energy; ++idx)
        {
            // Set ranges
            flast_cosine_profile[idx]->GetXaxis()->SetRangeUser(lidx, 1);
            flast_cosine_err[idx]->GetXaxis()->SetRangeUser(lidx, 1);
            sumrms_cosine_profile[idx]->GetXaxis()->SetRangeUser(lidx, 1);
            sumrms_cosine_err[idx]->GetXaxis()->SetRangeUser(lidx, 1);
            
            // Profile canvas
            flast_full_profiles->cd(idx+1);
            gStyle->SetOptStat(0);
            flast_cosine_profile[idx]->Draw();
            flast_profile_fit_func[idx]->Draw("same");
            sumrms_full_profiles->cd(idx+1);
            gStyle->SetOptStat(0);
            sumrms_cosine_profile[idx]->Draw();
            sumrms_profile_fit_func[idx]->Draw("same");

            // Profile error canvas
            flast_full_profiles_err->cd(idx+1);
            gStyle->SetOptStat(0);
            flast_cosine_err[idx]->Draw();
            flast_profile_err_fit_func[idx]->Draw("same");
            sumrms_full_profiles_err->cd(idx+1);
            gStyle->SetOptStat(0);
            sumrms_cosine_err[idx]->Draw();
            sumrms_profile_err_fit_func[idx]->Draw("same");
        }

        flast_full_profiles->Write();
        flast_full_profiles_err->Write();
        sumrms_full_profiles->Write();
        sumrms_full_profiles_err->Write();
    }

    // Write fit summary
    if (verbose)
        std::cout << "\nWriting fit summary [" << fit_summary_output << "]\n";
    FullFitSummary(
        flast_profile_fit_func,
        flast_profile_err_fit_func,
        sumrms_profile_fit_func,
        sumrms_profile_err_fit_func,
        fit_summary_output,
        "fit_summary");

}


void NUDFit(
    const char* full_histo_path, 
    const char* output)
{
    TFile *input_file = TFile::Open(full_histo_path, "READ");
    if (input_file->IsZombie())
    {
        std::cerr << "\n\nError opening input file: [" << full_histo_path << "]\n\n";
        exit(100);
    }

    auto h_nud_max = static_cast<TH1D*>(input_file->Get("Preselection/NUD/h_NUD_ps_total_adc"));
    h_nud_max->SetDirectory(0);
    input_file->Close();

    TF1 fitfunc("fitfunc", "log10([0]+[1]*x+[2]*pow(x,2))", 200, 10000);
    h_nud_max->Fit("fitfunc", "LR");

    TFile *output_file = TFile::Open(output, "RECREATE");
    if (output_file->IsZombie())
    {
        std::cerr << "\n\nError writing output file: [" << output << "]\n\n";
        exit(100);
    }

    h_nud_max->Write();

    output_file->Close();

}