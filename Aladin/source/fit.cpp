#include "DAMPE_geo_structure.h"
#include "fit.h"
#include "utils.h"

#include <map>
#include <vector>
#include <memory>

#include "TF1.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <ROOT/RDataFrame.hxx>

inline void extract_layer_lambda(
    std::vector<int> &best_histos_idx,
    std::vector<double> &best_lambda, 
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> in_gaus_histos,
    const double lambda_start,
    const double lambda_step,
    const int lambda_num)
{
    auto skew_limit = 0.1;
    std::vector<double> goodness (DAMPE_bgo_nLayers, 999);
    for (int lambda_idx=0; lambda_idx<=lambda_num; ++lambda_idx)
    {
        for (unsigned int ly_idx=0; ly_idx<DAMPE_bgo_nLayers; ++ly_idx)
        {
            auto xmin = in_gaus_histos[lambda_idx][ly_idx]->GetXaxis()->GetXmin();
            auto xmax = in_gaus_histos[lambda_idx][ly_idx]->GetXaxis()->GetXmax();
            std::unique_ptr<TF1> fitfunc = std::make_unique<TF1>("fitfunc", "gaus", xmin, xmax);
            fitfunc->SetNpx(10000);
            in_gaus_histos[lambda_idx][ly_idx]->Fit("fitfunc", "QW");
            auto skew = fabs(in_gaus_histos[lambda_idx][ly_idx]->GetSkewness());
            auto chi2 = fitfunc->GetChisquare();
            auto dof = fitfunc->GetNDF();

            if ((in_gaus_histos[lambda_idx][ly_idx]->GetRMS()/skew)>0.01 && (in_gaus_histos[lambda_idx][ly_idx]->GetRMS()/skew)<100.0 && skew<skew_limit)
            {
                double tmp_goodness = chi2/dof;
                if (goodness[ly_idx]==999)
                {
                    best_lambda[ly_idx] = lambda_start + lambda_idx*lambda_step;
                    best_histos_idx[ly_idx] = lambda_idx;
                    goodness[ly_idx] = tmp_goodness;
                }
                else
                {
                    if (tmp_goodness<goodness[ly_idx])
                    {
                        best_lambda[ly_idx] = lambda_start + lambda_idx*lambda_step;
                        best_histos_idx[ly_idx] = lambda_idx;
                        goodness[ly_idx] = tmp_goodness;
                    }
                }
            }
        }
    }
}

inline void extract_lambda(
    int &best_histo_idx, 
    double &best_lambda, 
    std::vector<ROOT::RDF::RResultPtr<TH1D>> in_gaus_histos,
    const double lambda_start,
    const double lambda_step,
    const int lambda_num)
{
    double goodness = 999;
    auto skew_limit = 0.1;
    for (int lambda_idx=0; lambda_idx<=lambda_num; ++lambda_idx)
    {
        auto xmin = in_gaus_histos[lambda_idx]->GetXaxis()->GetXmin();
        auto xmax = in_gaus_histos[lambda_idx]->GetXaxis()->GetXmax();
        std::unique_ptr<TF1> fitfunc = std::make_unique<TF1>("fitfunc", "gaus", xmin, xmax);
        fitfunc->SetNpx(10000);
        in_gaus_histos[lambda_idx]->Fit("fitfunc", "QW");
        auto skew = fabs(in_gaus_histos[lambda_idx]->GetSkewness());
        auto chi2 = fitfunc->GetChisquare();
        auto dof = fitfunc->GetNDF();

        if ((in_gaus_histos[lambda_idx]->GetRMS()/skew)>0.01 && (in_gaus_histos[lambda_idx]->GetRMS()/skew)<100.0 && skew<skew_limit)
        {
            double tmp_goodness = chi2/dof;
            if (goodness==999)
            {
                best_lambda = lambda_start + lambda_idx*lambda_step;
                best_histo_idx = lambda_idx;
                goodness = tmp_goodness;
            }
            else
            {
                if (tmp_goodness<goodness)
                {
                    best_lambda = lambda_start + lambda_idx*lambda_step;
                    best_histo_idx = lambda_idx;
                    goodness = tmp_goodness;
                }
            }
        }    
    }
}

void fit(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    std::shared_ptr<lambda_config> _lambda_config,
    const double _entries,
    unsigned int focus_energybin,
    const std::string outputPath,
    const std::string regularize_tree_path,
    const bool verbose,
    const unsigned int threads,
    const bool _mc)
{
    ROOT::EnableImplicitMT(threads);
    ROOT::RDataFrame _data_fr(*evtch);

    // Extract the energy binning
    auto energy_binning = _energy_config->GetEnergyBinning();

    auto rms_lambda_values = _lambda_config->GetRMSLambdaStruct();
    auto sumrms_lambda_values = _lambda_config->GetSumRMSLambdaStruct();
    auto elf_lambda_values = _lambda_config->GetELFLambdaStruct();
    auto ell_lambda_values = _lambda_config->GetELLLambdaStruct();
    auto xtrl_lambda_values = _lambda_config->GetXTRLLambdaStruct();

    if (verbose)
    {   
        _config->PrintActiveFilters();
        _energy_config->PrintActiveFilters();
        std::cout << "Total number of events: " << _entries;
        _lambda_config->PrintLambdaSettings();
        std::cout << "\nAnalysis running..." << std::endl;
    }
    
    TFile* output_file = TFile::Open(outputPath.c_str(), "RECREATE");
    if (output_file->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outputPath << "]\n\n";
        exit(100);
    }

    if (verbose) std::cout << "\nBuilding RMS and LFS distributions..." << std::endl;
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>>  h_rmslayer_gauss (rms_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_sumrms_gauss (sumrms_lambda_values.num+1);
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>>  h_energyfrac_layer_gauss (elf_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_energyfrac_last_layer_gauss (ell_lambda_values.num+1);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_xtrl_gauss (xtrl_lambda_values.num+1);

    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_rmslayer_gauss_scale (rms_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_sumrms_gauss_scale (sumrms_lambda_values.num+1);
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_energyfrac_layer_gauss_scale (elf_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_energyfrac_last_layer_gauss_scale (ell_lambda_values.num+1);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_xtrl_gauss_scale (xtrl_lambda_values.num+1);

    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_rmslayer_gauss_norm (rms_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_sumrms_gauss_norm (sumrms_lambda_values.num+1);
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_energyfrac_layer_gauss_norm (elf_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_energyfrac_last_layer_gauss_norm (ell_lambda_values.num+1);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_xtrl_gauss_norm (xtrl_lambda_values.num+1);
    
    double lambda;
    auto bin_filter = [focus_energybin](int energy_bin) -> bool { return energy_bin == (int)focus_energybin; };

    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {   
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx; 
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {   
            auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
            h_rmslayer_gauss[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"rmslayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {   
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx;   
        auto map_filter = [lambda](std::map<double, double> map_gauss) -> double { return map_gauss[lambda]; };
        h_sumrms_gauss[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"sumrms_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
    }

    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {   
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx; 
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
            h_energyfrac_layer_gauss[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fraclayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=ell_lambda_values.num; ++l_idx)
    {   
        lambda = ell_lambda_values.start + ell_lambda_values.step*l_idx; 
        auto map_filter = [lambda](std::map<double, double> map_gauss) -> double { return map_gauss[lambda]; };
        h_energyfrac_last_layer_gauss[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fraclastlayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
    }

    for (int l_idx=0; l_idx<=xtrl_lambda_values.num; ++l_idx)
    {   
        lambda = xtrl_lambda_values.start + xtrl_lambda_values.step*l_idx; 
        auto map_filter = [lambda](std::map<double, double> map_gauss) -> double { return map_gauss[lambda]; };
        h_xtrl_gauss[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"xtrl_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
    }

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/RMS")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/RMS")).c_str());
    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto h_name = std::string("h_rms_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);            
            h_rmslayer_gauss[l_idx][ly]->SetName(h_name.c_str());
            h_rmslayer_gauss[l_idx][ly]->GetXaxis()->SetTitle("RMS_{#lambda}");
            h_rmslayer_gauss[l_idx][ly]->Write();
        }
    }
    
    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/SumRMS")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/SumRMS")).c_str());
    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto h_name = std::string("h_sumrms_lambda_") + str_lambda;            
        h_sumrms_gauss[l_idx]->SetName(h_name.c_str());
        h_sumrms_gauss[l_idx]->GetXaxis()->SetTitle("SumRMS_{#lambda}");
        h_sumrms_gauss[l_idx]->Write();
    }
    
    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELF")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELF")).c_str());
    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto h_name = std::string("h_fraclayer_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly); 
            h_energyfrac_layer_gauss[l_idx][ly]->SetName(h_name.c_str());
            h_energyfrac_layer_gauss[l_idx][ly]->GetXaxis()->SetTitle("ELF_{#lambda}");
            h_energyfrac_layer_gauss[l_idx][ly]->Write();
        }
    }

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELL")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELL")).c_str());
    for (int l_idx=0; l_idx<=ell_lambda_values.num; ++l_idx)
    {
        lambda = ell_lambda_values.start + ell_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto h_name = std::string("h_fraclastlayer_lambda_") + str_lambda; 
        h_energyfrac_last_layer_gauss[l_idx]->SetName(h_name.c_str());
        h_energyfrac_last_layer_gauss[l_idx]->GetXaxis()->SetTitle("ELL_{#lambda}");
        h_energyfrac_last_layer_gauss[l_idx]->Write();
    }

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/XTRL")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/XTRL")).c_str());
    for (int l_idx=0; l_idx<=xtrl_lambda_values.num; ++l_idx)
    {
        lambda = xtrl_lambda_values.start + xtrl_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto h_name = std::string("h_xtrl_lambda_") + str_lambda; 
        h_xtrl_gauss[l_idx]->SetName(h_name.c_str());
        h_xtrl_gauss[l_idx]->GetXaxis()->SetTitle("XTRL_{#lambda}");
        h_xtrl_gauss[l_idx]->Write();
    }

    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx;
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {   
            auto map_filter_scale = [lambda, l_idx, ly, &h_rmslayer_gauss](std::map<double, std::vector<double>> map_gauss) -> double 
            {
                auto new_val = map_gauss[lambda][ly];
                auto hsigma = h_rmslayer_gauss[l_idx][ly]->GetRMS();
                if (hsigma) new_val /= hsigma;
                return new_val;
            };
            h_rmslayer_gauss_scale[l_idx][ly]  = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                                .Define("mapval", map_filter_scale, {"rmslayer_gauss"})
                                                                .Histo1D<double, double>("mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx;
        auto map_filter_scale = [lambda, l_idx, &h_sumrms_gauss](std::map<double, double> map_gauss) -> double 
        {   
            auto new_val = map_gauss[lambda];
            auto hsigma = h_sumrms_gauss[l_idx]->GetRMS();
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        h_sumrms_gauss_scale[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("mapval", map_filter_scale, {"sumrms_gauss"})
                                                    .Histo1D<double, double>("mapval", "simu_energy_w_corr");
    }

    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx;
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {   
            auto map_filter_scale = [lambda, l_idx, ly, &h_energyfrac_layer_gauss](std::map<double, std::vector<double>> map_gauss) -> double 
            {   
                auto new_val = map_gauss[lambda][ly];
                auto hsigma = h_energyfrac_layer_gauss[l_idx][ly]->GetRMS();
                if (hsigma) new_val /= hsigma;
                return new_val;
            };
            h_energyfrac_layer_gauss_scale[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                        .Define("mapval", map_filter_scale, {"fraclayer_gauss"})
                                                        .Histo1D<double, double>("mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=ell_lambda_values.num; ++l_idx)
    {
        lambda = ell_lambda_values.start + ell_lambda_values.step*l_idx;
        auto map_filter_scale = [lambda, l_idx, &h_energyfrac_last_layer_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hsigma = h_energyfrac_last_layer_gauss[l_idx]->GetRMS();
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        h_energyfrac_last_layer_gauss_scale[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("mapval", map_filter_scale, {"fraclastlayer_gauss"})
                                                    .Histo1D<double, double>("mapval", "simu_energy_w_corr");
    }

    for (int l_idx=0; l_idx<=xtrl_lambda_values.num; ++l_idx)
    {
        lambda = xtrl_lambda_values.start + xtrl_lambda_values.step*l_idx;
        auto map_filter_scale = [lambda, l_idx, &h_xtrl_gauss](std::map<double, double> map_gauss) -> double 
        {   
            auto new_val = map_gauss[lambda];
            auto hsigma = h_xtrl_gauss[l_idx]->GetRMS();
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        h_xtrl_gauss_scale[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                            .Define("mapval", map_filter_scale, {"xtrl_gauss"})
                                            .Histo1D<double, double>("mapval", "simu_energy_w_corr");
    }

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/RMS_scale")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/RMS_scale")).c_str());
    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto h_name = std::string("h_rms_scale_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);            
            h_rmslayer_gauss_scale[l_idx][ly]->SetName(h_name.c_str());
            h_rmslayer_gauss_scale[l_idx][ly]->GetXaxis()->SetTitle("RMS_{#lambda}");
            h_rmslayer_gauss_scale[l_idx][ly]->Write();
        }
    }
    
    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/SumRMS_scale")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/SumRMS_scale")).c_str());
    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto h_name = std::string("h_sumrms_scale_lambda_") + str_lambda;            
        h_sumrms_gauss_scale[l_idx]->SetName(h_name.c_str());
        h_sumrms_gauss_scale[l_idx]->GetXaxis()->SetTitle("SumRMS_{#lambda}");
        h_sumrms_gauss_scale[l_idx]->Write();
    }
    
    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELF_scale")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELF_scale")).c_str());
    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto h_name = std::string("h_fraclayer_scale_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly); 
            h_energyfrac_layer_gauss_scale[l_idx][ly]->SetName(h_name.c_str());
            h_energyfrac_layer_gauss_scale[l_idx][ly]->GetXaxis()->SetTitle("ELF_{#lambda}");
            h_energyfrac_layer_gauss_scale[l_idx][ly]->Write();
        }
    }

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELL_scale")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELL_scale")).c_str());
    for (int l_idx=0; l_idx<=ell_lambda_values.num; ++l_idx)
    {
        lambda = ell_lambda_values.start + ell_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto h_name = std::string("h_fraclastlayer_scale_lambda_") + str_lambda; 
        h_energyfrac_last_layer_gauss_scale[l_idx]->SetName(h_name.c_str());
        h_energyfrac_last_layer_gauss_scale[l_idx]->GetXaxis()->SetTitle("ELL_{#lambda}");
        h_energyfrac_last_layer_gauss_scale[l_idx]->Write();
    }

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/XTRL_scale")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/XTRL_scale")).c_str());
    for (int l_idx=0; l_idx<=xtrl_lambda_values.num; ++l_idx)
    {
        lambda = xtrl_lambda_values.start + xtrl_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto h_name = std::string("h_xtrl_scale_lambda_") + str_lambda; 
        h_xtrl_gauss_scale[l_idx]->SetName(h_name.c_str());
        h_xtrl_gauss_scale[l_idx]->GetXaxis()->SetTitle("XTRL_{#lambda}");
        h_xtrl_gauss_scale[l_idx]->Write();
    }

    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {   
            auto map_filter_scale = [lambda, l_idx, ly, &h_rmslayer_gauss](std::map<double, std::vector<double>> map_gauss) -> double 
            {
                auto new_val = map_gauss[lambda][ly];
                auto hsigma = h_rmslayer_gauss[l_idx][ly]->GetRMS();
                if (hsigma) new_val /= hsigma;
                return new_val;
            };
            auto map_filter_shift = [l_idx, ly, &h_rmslayer_gauss_scale](const double map_gauss) -> double 
            { 
                auto new_val = map_gauss;
                auto hmean = h_rmslayer_gauss_scale[l_idx][ly]->GetMean();
                if (hmean>0) new_val -= hmean;
                else new_val += abs(hmean);
                return new_val;
            };
            auto h_name = std::string("h_rms_norm_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
            auto h_title = std::string("Normalized RMS - #lambda ") + str_lambda + std::string(" - layer ") + std::to_string(ly);
            h_rmslayer_gauss_norm[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                        .Define("mapval_scale", map_filter_scale, {"rmslayer_gauss"})
                                                        .Define("mapval_shift", map_filter_shift, {"mapval_scale"})
                                                        .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval_shift", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto map_filter_scale = [lambda, l_idx, &h_sumrms_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hsigma = h_sumrms_gauss[l_idx]->GetRMS();
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto map_filter_shift = [l_idx, &h_sumrms_gauss_scale](double map_gauss) -> double 
        { 
            auto new_val = map_gauss;
            auto hmean = h_sumrms_gauss_scale[l_idx]->GetMean();
            if (hmean>0) new_val -= hmean;
            else new_val += abs(hmean);
            return new_val;
        };
        auto h_name = std::string("h_sumrms_norm_lambda_") + str_lambda;
        auto h_title = std::string("Normalized SumRMS - #lambda ") + str_lambda;
        h_sumrms_gauss_norm[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("mapval_scale", map_filter_scale, {"sumrms_gauss"})
                                                    .Define("mapval_shift", map_filter_shift, {"mapval_scale"})
                                                    .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval_shift", "simu_energy_w_corr");
    }

    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {   
            auto map_filter_scale = [lambda, l_idx, ly, &h_energyfrac_layer_gauss](std::map<double, std::vector<double>> map_gauss) -> double 
            {   
                auto new_val = map_gauss[lambda][ly];
                auto hsigma = h_energyfrac_layer_gauss[l_idx][ly]->GetRMS();
                if (hsigma) new_val /= hsigma;
                return new_val;
            };
            auto map_filter_shift = [l_idx, ly, &h_energyfrac_layer_gauss_scale](const double map_gauss) -> double 
            { 
                auto new_val = map_gauss;
                auto hmean = h_energyfrac_layer_gauss_scale[l_idx][ly]->GetMean();
                if (hmean>0) new_val -= hmean;
                else new_val += abs(hmean);
                return new_val;
            };
            auto h_name = std::string("h_fraclayer_norm_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
            auto h_title = std::string("Normalized ELF - #lambda ") + str_lambda + std::string(" - layer ") + std::to_string(ly);
            h_energyfrac_layer_gauss_norm[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                        .Define("mapval_scale", map_filter_scale, {"fraclayer_gauss"})
                                                        .Define("mapval_shift", map_filter_shift, {"mapval_scale"})
                                                        .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval_shift", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=ell_lambda_values.num; ++l_idx)
    {
        lambda = ell_lambda_values.start + ell_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto map_filter_scale = [lambda, l_idx, &h_energyfrac_last_layer_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hsigma = h_energyfrac_last_layer_gauss[l_idx]->GetRMS();
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto map_filter_shift = [l_idx, &h_energyfrac_last_layer_gauss_scale](const double map_gauss) -> double 
        { 
            auto new_val = map_gauss;
            auto hmean = h_energyfrac_last_layer_gauss_scale[l_idx]->GetMean();
            if (hmean>0) new_val -= hmean;
            else new_val += abs(hmean);
            return new_val;
        };
        auto h_name = std::string("h_fraclastlayer_norm_lambda_") + str_lambda;
        auto h_title = std::string("Normalized ELL - #lambda ") + str_lambda;
        h_energyfrac_last_layer_gauss_norm[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("mapval_scale", map_filter_scale, {"fraclastlayer_gauss"})
                                                    .Define("mapval_shift", map_filter_shift, {"mapval_scale"})
                                                    .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval_shift", "simu_energy_w_corr");
    }
    
    for (int l_idx=0; l_idx<=xtrl_lambda_values.num; ++l_idx)
    {
        lambda = xtrl_lambda_values.start + xtrl_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto map_filter_scale = [lambda, l_idx, &h_xtrl_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hsigma = h_xtrl_gauss[l_idx]->GetRMS();
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto map_filter_shift = [l_idx, &h_xtrl_gauss_scale](const double map_gauss) -> double 
        { 
            auto new_val = map_gauss;
            auto hmean = h_xtrl_gauss_scale[l_idx]->GetMean();
            if (hmean>0) new_val -= hmean;
            else new_val += abs(hmean);
            return new_val;
        };
        auto h_name = std::string("h_xtrl_norm_lambda_") + str_lambda;
        auto h_title = std::string("Normalized XTRL - #lambda ") + str_lambda;
        h_xtrl_gauss_norm[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                            .Define("mapval_scale", map_filter_scale, {"xtrl_gauss"})
                                            .Define("mapval_shift", map_filter_shift, {"mapval_scale"})
                                            .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval_shift", "simu_energy_w_corr");
    }

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/RMS_norm")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/RMS_norm")).c_str());
    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            h_rmslayer_gauss_norm[l_idx][ly]->GetXaxis()->SetTitle("RMS_{#lambda}");
            h_rmslayer_gauss_norm[l_idx][ly]->Write();
        }
    
    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/SumRMS_norm")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/SumRMS_norm")).c_str());
    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {
        h_sumrms_gauss_norm[l_idx]->GetXaxis()->SetTitle("SumRMS_{#lambda}");
        h_sumrms_gauss_norm[l_idx]->Write();
    }

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELF_norm")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELF_norm")).c_str());
    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            h_energyfrac_layer_gauss_norm[l_idx][ly]->GetXaxis()->SetTitle("ELF_{#lambda}");
            h_energyfrac_layer_gauss_norm[l_idx][ly]->Write();
        }

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELL_norm")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELL_norm")).c_str());
    for (int l_idx=0; l_idx<=ell_lambda_values.num; ++l_idx)
    {
        h_energyfrac_last_layer_gauss_norm[l_idx]->GetXaxis()->SetTitle("ELL_{#lambda}");
        h_energyfrac_last_layer_gauss_norm[l_idx]->Write();
    }
    
    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/XTRL_norm")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/XTRL_norm")).c_str());
    for (int l_idx=0; l_idx<=xtrl_lambda_values.num; ++l_idx)
    {
        h_xtrl_gauss_norm[l_idx]->GetXaxis()->SetTitle("XTRL_{#lambda}");
        h_xtrl_gauss_norm[l_idx]->Write();
    }

    // Find best lambda values
    if (verbose) std::cout << "\nFinding best lambda values..." << std::endl;
    std::vector<int> best_rms_hist_idx (DAMPE_bgo_nLayers, 0);
    int best_sumrms_hist_idx = 0;
    std::vector<int> best_fraclayer_hist_idx (DAMPE_bgo_nLayers, 0);
    int best_fraclast_hist_idx = 0;
    int best_xtrl_hist_idx = 0;
    
    auto lambda_select = best_lambda();
    lambda_select.rms = std::vector<double> (DAMPE_bgo_nLayers, rms_lambda_values.start);
    lambda_select.sumrms = sumrms_lambda_values.start;
    lambda_select.fraclayer = std::vector<double> (DAMPE_bgo_nLayers, elf_lambda_values.start);
    lambda_select.fraclast = ell_lambda_values.start;
    lambda_select.xtrl = xtrl_lambda_values.start;

    extract_layer_lambda(
        best_rms_hist_idx, 
        lambda_select.rms, 
        h_rmslayer_gauss_norm, 
        rms_lambda_values.start, 
        rms_lambda_values.step, 
        rms_lambda_values.num);
    
    extract_layer_lambda(
        best_fraclayer_hist_idx, 
        lambda_select.fraclayer, 
        h_energyfrac_layer_gauss_norm, 
        elf_lambda_values.start, 
        elf_lambda_values.step, 
        elf_lambda_values.num);
    
    extract_lambda(
        best_sumrms_hist_idx,
        lambda_select.sumrms,
        h_sumrms_gauss_norm,
        sumrms_lambda_values.start,
        sumrms_lambda_values.step,
        sumrms_lambda_values.num);

    extract_lambda(
        best_fraclast_hist_idx,
        lambda_select.fraclast,
        h_energyfrac_last_layer_gauss_norm,
        ell_lambda_values.start,
        ell_lambda_values.step,
        ell_lambda_values.num);
    
    extract_lambda(
        best_xtrl_hist_idx,
        lambda_select.xtrl,
        h_xtrl_gauss_norm,
        xtrl_lambda_values.start,
        xtrl_lambda_values.step,
        xtrl_lambda_values.num);
    
    std::unique_ptr<TCanvas> rms_bestfit = std::make_unique<TCanvas>("rms_bestfit", "RMS bestfit");
    std::unique_ptr<TCanvas> sumrms_bestfit = std::make_unique<TCanvas>("sumrms_bestfit", "SumRMS bestfit");
    std::unique_ptr<TCanvas> elf_bestfit = std::make_unique<TCanvas>("elf_bestfit", "ELF bestfit");
    std::unique_ptr<TCanvas> ell_bestfit = std::make_unique<TCanvas>("ell_bestfit", "ELL bestfit");
    std::unique_ptr<TCanvas> xtrl_bestfit = std::make_unique<TCanvas>("xtrl_bestfit", "XTRL bestfit");

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/Canvas")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/Canvas")).c_str());
    
    rms_bestfit->Divide(7,2);
    for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
    {
        rms_bestfit->cd(ly+1);
        h_rmslayer_gauss_norm[best_rms_hist_idx[ly]][ly]->Draw();
    }
    rms_bestfit->cd(0);
    rms_bestfit->Write();

    sumrms_bestfit->cd();
    h_sumrms_gauss_norm[best_sumrms_hist_idx]->Draw();
    sumrms_bestfit->Write();

    elf_bestfit->Divide(7,2);
    for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
    {
        elf_bestfit->cd(ly+1);
        h_energyfrac_layer_gauss_norm[best_fraclayer_hist_idx[ly]][ly]->Draw();
    }
    elf_bestfit->cd(0);
    elf_bestfit->Write();

    ell_bestfit->cd();
    h_energyfrac_last_layer_gauss_norm[best_fraclast_hist_idx]->Draw();
    ell_bestfit->Write();

    xtrl_bestfit->cd();
    h_xtrl_gauss_norm[best_xtrl_hist_idx]->Draw();
    xtrl_bestfit->Write();

    output_file->Close();

    // Build summary TTree
    std::string tree_vars_file = outputPath.substr(0, outputPath.find_last_of("/")) + std::string("/lambda_corrections.root");
    TFile* corrections_file = TFile::Open(tree_vars_file.c_str(), "RECREATE");
    if (corrections_file->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << tree_vars_file << "]\n\n";
        exit(100);
    }

    TTree corrections_tree("corrections_tree", "Lambda corrections");

    corrections_tree.Branch("energy_bin", &focus_energybin, "energy_bin/i");
    corrections_tree.Branch("best_rms_lambda", &lambda_select.rms);
    corrections_tree.Branch("best_sumrms_lambda", &lambda_select.sumrms, "best_sumrms_lambda/D");
    corrections_tree.Branch("best_fraclayer_lambda", &lambda_select.fraclayer);
    corrections_tree.Branch("best_fraclast_lambda", &lambda_select.fraclast, "best_fraclast_lambda/D");
    corrections_tree.Branch("best_xtrl_lambda", &lambda_select.xtrl, "best_xtrl_lambda/D");

    // Build vectors for Mean and RMS used to obtain normalized distributions
    std::vector<double> rms_norm_mean (DAMPE_bgo_nLayers, 0);
    std::vector<double> rms_norm_rms (DAMPE_bgo_nLayers, 0);
    double sumrms_norm_mean = h_sumrms_gauss_scale[best_sumrms_hist_idx]->GetMean();
    double sumrms_norm_rms = h_sumrms_gauss[best_sumrms_hist_idx]->GetRMS();
    std::vector<double> fraclayer_norm_mean (DAMPE_bgo_nLayers, 0);
    std::vector<double> fraclayer_norm_rms (DAMPE_bgo_nLayers, 0);
    double fraclast_norm_mean = h_energyfrac_last_layer_gauss_scale[best_fraclast_hist_idx]->GetMean();
    double fraclast_norm_rms = h_energyfrac_last_layer_gauss[best_fraclast_hist_idx]->GetRMS();
    double xtrl_norm_mean = h_xtrl_gauss_scale[best_xtrl_hist_idx]->GetMean();
    double xtrl_norm_rms = h_xtrl_gauss[best_xtrl_hist_idx]->GetRMS();
    
    for (int l_idx=0; l_idx<DAMPE_bgo_nLayers; ++l_idx)
    {
        rms_norm_mean[l_idx] = h_rmslayer_gauss_scale[best_rms_hist_idx[l_idx]][l_idx]->GetMean();
        rms_norm_rms[l_idx] = h_rmslayer_gauss[best_rms_hist_idx[l_idx]][l_idx]->GetRMS();
        fraclayer_norm_mean[l_idx] = h_energyfrac_layer_gauss[best_fraclayer_hist_idx[l_idx]][l_idx]->GetRMS();
        fraclayer_norm_rms[l_idx] = h_energyfrac_layer_gauss_scale[best_fraclayer_hist_idx[l_idx]][l_idx]->GetMean();
    }
    
    corrections_tree.Branch("rms_norm_mean", &rms_norm_mean);
    corrections_tree.Branch("rms_norm_rms", &rms_norm_rms);
    corrections_tree.Branch("sumrms_norm_mean", &sumrms_norm_mean, "sumrms_norm_mean/D");
    corrections_tree.Branch("sumrms_norm_rms", &sumrms_norm_rms, "sumrms_norm_rms/D");
    corrections_tree.Branch("fraclayer_norm_mean", &fraclayer_norm_mean);
    corrections_tree.Branch("fraclayer_norm_rms", &fraclayer_norm_rms);
    corrections_tree.Branch("fraclast_norm_mean", &fraclast_norm_mean, "fraclast_norm_mean/D");
    corrections_tree.Branch("fraclast_norm_rms", &fraclast_norm_rms, "fraclast_norm_rms/D");
    corrections_tree.Branch("xtrl_norm_mean", &xtrl_norm_mean, "xtrl_norm_mean/D");
    corrections_tree.Branch("xtrl_norm_rms", &xtrl_norm_rms, "xtrl_norm_rms/D");

    corrections_tree.Fill();
    corrections_tree.Write();
    corrections_file->Close();

    if (verbose)
    {
        std::cout << "\n\nOutput file has been written... [" << outputPath << "]\n";
        std::cout << "Output corrections file has been written... [" << tree_vars_file << "]\n";
    }

}