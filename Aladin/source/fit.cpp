#include "DAMPE_geo_structure.h"
#include "fit.h"

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
    std::vector<ROOT::RDF::RResultPtr<TH1D>> &best_histos, 
    std::vector<double> &best_lambda, 
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> in_gaus_histos,
    const double lambda_start,
    const double lambda_step,
    const int lambda_num)
{
    int min_stats = 5;
    std::vector<double> goodness (DAMPE_bgo_nLayers, 999);
    for (int lambda_idx=0; lambda_idx<=lambda_num; ++lambda_idx)
    {
        for (unsigned int ly_idx=0; ly_idx<DAMPE_bgo_nLayers; ++ly_idx)
        {
            if (in_gaus_histos[lambda_idx][ly_idx]->GetEntries()>min_stats)
            {
                auto xmin = in_gaus_histos[lambda_idx][ly_idx]->GetXaxis()->GetXmin();
                auto xmax = in_gaus_histos[lambda_idx][ly_idx]->GetXaxis()->GetXmax();
                std::unique_ptr<TF1> fitfunc = std::make_unique<TF1>("fitfunc", "gaus", xmin, xmax);
                fitfunc->SetNpx(10000);
                in_gaus_histos[lambda_idx][ly_idx]->Fit("fitfunc", "QW");
                auto chi2 = fitfunc->GetChisquare();
                auto dof = fitfunc->GetNDF();

                if (dof)
                {
                    double tmp_goodness = chi2/dof;
                    if (goodness[ly_idx]==999)
                    {
                        best_lambda[ly_idx] = lambda_start + lambda_idx*lambda_step;
                        best_histos[ly_idx] = in_gaus_histos[lambda_idx][ly_idx];
                        goodness[ly_idx] = tmp_goodness;
                    }
                    else
                    {
                        if (abs(tmp_goodness-1) < abs(goodness[ly_idx]-1))
                        {
                            best_lambda[ly_idx] = lambda_start + lambda_idx*lambda_step;
                            best_histos[ly_idx] = in_gaus_histos[lambda_idx][ly_idx];
                            goodness[ly_idx] = tmp_goodness;
                        }
                    }
                }

            }
        }
    }
    
}

inline void extract_lambda(
    ROOT::RDF::RResultPtr<TH1D> &best_histo, 
    double &best_lambda, 
    std::vector<ROOT::RDF::RResultPtr<TH1D>> in_gaus_histos,
    const double lambda_start,
    const double lambda_step,
    const int lambda_num)
{
    int min_stats = 5;
    double goodness = 999;
    for (int lambda_idx=0; lambda_idx<=lambda_num; ++lambda_idx)
    {
        if (in_gaus_histos[lambda_idx]->GetEntries()>min_stats)
        {
            auto xmin = in_gaus_histos[lambda_idx]->GetXaxis()->GetXmin();
            auto xmax = in_gaus_histos[lambda_idx]->GetXaxis()->GetXmax();
            std::unique_ptr<TF1> fitfunc = std::make_unique<TF1>("fitfunc", "gaus", xmin, xmax);
            fitfunc->SetNpx(10000);
            in_gaus_histos[lambda_idx]->Fit("fitfunc", "QW");
            auto chi2 = fitfunc->GetChisquare();
            auto dof = fitfunc->GetNDF();

            if (dof)
            {
                double tmp_goodness = chi2/dof;
                if (goodness==999)
                {
                    best_lambda = lambda_start + lambda_idx*lambda_step;
                    best_histo = in_gaus_histos[lambda_idx];
                    goodness = tmp_goodness;
                }
                else
                {
                    if (abs(tmp_goodness-1) < abs(goodness-1))
                    {
                        best_lambda = lambda_start + lambda_idx*lambda_step;
                        best_histo = in_gaus_histos[lambda_idx];
                        goodness = tmp_goodness;
                    }
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
    auto elf_ang_lambda_values = _lambda_config->GetELFAngLambdaStruct();

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
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_sumrmslayer_gauss (sumrms_lambda_values.num+1);
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>>  h_energyfrac_layer_gauss (elf_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_energyfrac_last_layer_gauss (elf_ang_lambda_values.num+1);

    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_rmslayer_gauss_norm (rms_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_sumrmslayer_gauss_norm (sumrms_lambda_values.num+1);
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_energyfrac_layer_gauss_norm (elf_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_energyfrac_last_layer_gauss_norm (elf_ang_lambda_values.num+1);
    
    double lambda;
    auto bin_filter = [focus_energybin](int energy_bin) -> bool { return energy_bin == (int)focus_energybin; };

    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {   
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx; 
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {   
            auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
            h_rmslayer_gauss[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"rmsLayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {   
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx;   
        auto map_filter = [lambda](std::map<double, double> map_gauss) -> double { return map_gauss[lambda]; };
        h_sumrmslayer_gauss[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"sumrmsLayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
    }

    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {   
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx; 
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
            h_energyfrac_layer_gauss[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fracLayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=elf_ang_lambda_values.num; ++l_idx)
    {   
        lambda = elf_ang_lambda_values.start + elf_ang_lambda_values.step*l_idx; 
        auto map_filter = [lambda](std::map<double, double> map_gauss) -> double { return map_gauss[lambda]; };
        h_energyfrac_last_layer_gauss[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fracLayer_ang_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
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
        h_sumrmslayer_gauss[l_idx]->SetName(h_name.c_str());
        h_sumrmslayer_gauss[l_idx]->GetXaxis()->SetTitle("SumRMS_{#lambda}");
        h_sumrmslayer_gauss[l_idx]->Write();
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
    for (int l_idx=0; l_idx<=elf_ang_lambda_values.num; ++l_idx)
    {
        lambda = elf_ang_lambda_values.start + elf_ang_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto h_name = std::string("h_fraclayer_ang_lambda_") + str_lambda; 
        h_energyfrac_last_layer_gauss[l_idx]->SetName(h_name.c_str());
        h_energyfrac_last_layer_gauss[l_idx]->GetXaxis()->SetTitle("ELL_{#lambda}");
        h_energyfrac_last_layer_gauss[l_idx]->Write();
    }

    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {   
            auto map_filter = [lambda, l_idx, ly, &h_rmslayer_gauss](std::map<double, std::vector<double>> map_gauss) -> double 
            { 
                auto new_val = map_gauss[lambda][ly];
                auto hmean = h_rmslayer_gauss[l_idx][ly]->GetMean();
                auto hsigma = h_rmslayer_gauss[l_idx][ly]->GetRMS();
                new_val -= hmean>0 ? hmean : -hmean;
                if (hsigma) new_val /= hsigma;
                return new_val;
            };
            auto h_name = std::string("h_rms_norm_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
            auto h_title = std::string("Normalized RMS - #lambda ") + str_lambda + std::string(" - layer ") + std::to_string(ly);
            h_rmslayer_gauss_norm[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                        .Define("mapval", map_filter, {"rmsLayer_gauss"})
                                                        .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto map_filter = [lambda, l_idx, &h_sumrmslayer_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hmean = h_sumrmslayer_gauss[l_idx]->GetMean();
            auto hsigma = h_sumrmslayer_gauss[l_idx]->GetRMS();
            new_val -= hmean>0 ? hmean : -hmean;
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto h_name = std::string("h_sumrms_norm_lambda_") + str_lambda;
        auto h_title = std::string("Normalized SumRMS - #lambda ") + str_lambda;
        h_sumrmslayer_gauss_norm[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("mapval", map_filter, {"sumrmsLayer_gauss"})
                                                    .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
    }

    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {   
            auto map_filter = [lambda, l_idx, ly, &h_energyfrac_layer_gauss](std::map<double, std::vector<double>> map_gauss) -> double 
            { 
                auto new_val = map_gauss[lambda][ly];
                auto hmean = h_energyfrac_layer_gauss[l_idx][ly]->GetMean();
                auto hsigma = h_energyfrac_layer_gauss[l_idx][ly]->GetRMS();
                new_val -= hmean>0 ? hmean : -hmean;
                if (hsigma) new_val /= hsigma;
                return new_val;
            };
            auto h_name = std::string("h_fraclayer_norm_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
            auto h_title = std::string("Normalized ELF - #lambda ") + str_lambda + std::string(" - layer ") + std::to_string(ly);
            h_energyfrac_layer_gauss_norm[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                        .Define("mapval", map_filter, {"fracLayer_gauss"})
                                                        .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=elf_ang_lambda_values.num; ++l_idx)
    {
        lambda = elf_ang_lambda_values.start + elf_ang_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto map_filter = [lambda, l_idx, &h_energyfrac_last_layer_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hmean = h_energyfrac_last_layer_gauss[l_idx]->GetMean();
            auto hsigma = h_energyfrac_last_layer_gauss[l_idx]->GetRMS();
            new_val -= hmean>0 ? hmean : -hmean;
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto h_name = std::string("h_fraclayer_ang_norm_lambda_") + str_lambda;
        auto h_title = std::string("Normalized ELFang - #lambda ") + str_lambda;
        h_energyfrac_last_layer_gauss_norm[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("mapval", map_filter, {"fracLayer_ang_gauss"})
                                                    .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
    }
    
    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/RMS_norm")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/RMS_norm")).c_str());
    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            h_rmslayer_gauss_norm[l_idx][ly]->Write();
            h_rmslayer_gauss_norm[l_idx][ly]->GetXaxis()->SetTitle("RMS_{#lambda}");
        }
    
    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/SumRMS_norm")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/SumRMS_norm")).c_str());
    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {
        h_sumrmslayer_gauss_norm[l_idx]->Write();
        h_sumrmslayer_gauss_norm[l_idx]->GetXaxis()->SetTitle("SumRMS_{#lambda}");
    }

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELF_norm")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELF_norm")).c_str());
    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            h_energyfrac_layer_gauss_norm[l_idx][ly]->Write();
            h_energyfrac_layer_gauss_norm[l_idx][ly]->GetXaxis()->SetTitle("ELF_{#lambda}");
        }

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELL_norm")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/ELL_norm")).c_str());
    for (int l_idx=0; l_idx<=elf_ang_lambda_values.num; ++l_idx)
    {
        h_energyfrac_last_layer_gauss_norm[l_idx]->Write();
        h_energyfrac_last_layer_gauss_norm[l_idx]->GetXaxis()->SetTitle("ELL_{#lambda}");
    }

    // Find best lambda values
    std::vector<ROOT::RDF::RResultPtr<TH1D>> best_rms_hist (DAMPE_bgo_nLayers);
    ROOT::RDF::RResultPtr<TH1D> best_sumrmshist;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> best_fraclayer_hist (DAMPE_bgo_nLayers);
    ROOT::RDF::RResultPtr<TH1D> best_fraclast_ang_hist;

    std::vector<double> best_rms_lambda (DAMPE_bgo_nLayers, 999);
    double best_sumrms_lambda;
    std::vector<double> best_fraclayer_lambda (DAMPE_bgo_nLayers, 999);
    double best_fraclast_lambda;
    
    extract_layer_lambda(
        best_rms_hist, 
        best_rms_lambda, 
        h_rmslayer_gauss_norm, 
        rms_lambda_values.start, 
        rms_lambda_values.step, 
        rms_lambda_values.num);
    
    extract_layer_lambda(
        best_fraclayer_hist, 
        best_fraclayer_lambda, 
        h_energyfrac_layer_gauss_norm, 
        elf_lambda_values.start, 
        elf_lambda_values.step, 
        elf_lambda_values.num);

    extract_lambda(
        best_sumrmshist,
        best_sumrms_lambda,
        h_sumrmslayer_gauss_norm,
        sumrms_lambda_values.start,
        sumrms_lambda_values.step,
        sumrms_lambda_values.num);

    extract_lambda(
        best_fraclast_ang_hist,
        best_fraclast_lambda,
        h_energyfrac_last_layer_gauss_norm,
        elf_ang_lambda_values.start,
        elf_ang_lambda_values.step,
        elf_ang_lambda_values.num);

    std::unique_ptr<TCanvas> rms_bestfit = std::make_unique<TCanvas>("rms_bestfit", "RMS bestfit");
    std::unique_ptr<TCanvas> sumrms_bestfit = std::make_unique<TCanvas>("sumrms_bestfit", "SumRMS bestfit");
    std::unique_ptr<TCanvas> elf_bestfit = std::make_unique<TCanvas>("elf_bestfit", "ELF bestfit");
    std::unique_ptr<TCanvas> ell_bestfit = std::make_unique<TCanvas>("ell_bestfit", "ELL bestfit");


    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/Canvas")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/Canvas")).c_str());
    
    rms_bestfit->Divide(7,2);
    for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
    {
        rms_bestfit->cd(ly+1);
        best_rms_hist[ly]->Draw();
    }
    rms_bestfit->cd(0);
    rms_bestfit->Write();

    sumrms_bestfit->cd();
    best_sumrmshist->Write();

    elf_bestfit->Divide(7,2);
    for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
    {
        elf_bestfit->cd(ly+1);
        best_fraclayer_hist[ly]->Draw();
    }
    elf_bestfit->cd(0);
    elf_bestfit->Write();

    ell_bestfit->cd();
    best_fraclast_ang_hist->Write();

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
    corrections_tree.Branch("best_rms_lambda", &best_rms_lambda);
    corrections_tree.Branch("best_sumrms_lambda", &best_sumrms_lambda, "best_sumrms_lambda/D");
    corrections_tree.Branch("best_fraclayer_lambda", &best_fraclayer_lambda);
    corrections_tree.Branch("best_fraclast_lambda", &best_fraclast_lambda, "best_fraclast_lambda/D");

    corrections_tree.Fill();
    corrections_tree.Write();
    corrections_file->Close();

    if (verbose)
    {
        std::cout << "\n\nOutput file has been written... [" << outputPath << "]\n";
        std::cout << "Output corrections file has been written... [" << tree_vars_file << "]\n";
    }
}