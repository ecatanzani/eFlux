#include "DAMPE_geo_structure.h"
#include "loglikelihood.h"

#include <map>
#include <vector>
#include <memory>

#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include <ROOT/RDataFrame.hxx>

void buildLogLikelihoodProfile(
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
        h_energyfrac_last_layer_gauss[l_idx]->SetName(h_name.c_str());
        h_energyfrac_last_layer_gauss[l_idx]->GetXaxis()->SetTitle("XTRL_{#lambda}");
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
                                                        .Define("mapval", map_filter, {"rmslayer_gauss"})
                                                        .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto map_filter = [lambda, l_idx, &h_sumrms_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hmean = h_sumrms_gauss[l_idx]->GetMean();
            auto hsigma = h_sumrms_gauss[l_idx]->GetRMS();
            new_val -= hmean>0 ? hmean : -hmean;
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto h_name = std::string("h_sumrms_norm_lambda_") + str_lambda;
        auto h_title = std::string("Normalized SumRMS - #lambda ") + str_lambda;
        h_sumrms_gauss_norm[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("mapval", map_filter, {"sumrms_gauss"})
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
                                                        .Define("mapval", map_filter, {"fraclayer_gauss"})
                                                        .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=ell_lambda_values.num; ++l_idx)
    {
        lambda = ell_lambda_values.start + ell_lambda_values.step*l_idx;
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
        auto h_name = std::string("h_fraclastlayer_norm_lambda_") + str_lambda;
        auto h_title = std::string("Normalized ELL - #lambda ") + str_lambda;
        h_energyfrac_last_layer_gauss_norm[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                            .Define("mapval", map_filter, {"fraclastlayer_gauss"})
                                                            .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
    }

    for (int l_idx=0; l_idx<=xtrl_lambda_values.num; ++l_idx)
    {
        lambda = xtrl_lambda_values.start + xtrl_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto map_filter = [lambda, l_idx, &h_xtrl_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hmean = h_xtrl_gauss[l_idx]->GetMean();
            auto hsigma = h_xtrl_gauss[l_idx]->GetRMS();
            new_val -= hmean>0 ? hmean : -hmean;
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto h_name = std::string("h_xtrl_norm_lambda_") + str_lambda;
        auto h_title = std::string("Normalized XTRL - #lambda ") + str_lambda;
        h_xtrl_gauss_norm[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                            .Define("mapval", map_filter, {"xtrl_gauss"})
                                            .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
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

    if (verbose) std::cout << "\nBuilding RMS and ELF likelihood profiles...\n";
    
    std::vector<std::vector<double>> likeprofile_rms (DAMPE_bgo_nLayers, std::vector<double> (rms_lambda_values.num+1));
    std::vector<double> likeprofile_sumrms (sumrms_lambda_values.num+1);
    std::vector<std::vector<double>> likeprofile_elf (DAMPE_bgo_nLayers, std::vector<double> (elf_lambda_values.num+1));
    std::vector<double> likeprofile_ell (ell_lambda_values.num+1);
    std::vector<double> likeprofile_xtrl (xtrl_lambda_values.num+1);
    
    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {   
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx; 
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
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
            auto loglike = [l_idx, ly, &h_rmslayer_gauss_norm](double value) -> double 
            {
                double pdf = TMath::Gaus(value, h_rmslayer_gauss_norm[l_idx][ly]->GetMean(), h_rmslayer_gauss_norm[l_idx][ly]->GetRMS(), true);
                double loolikelihood = pdf>0 ? std::log(pdf) : 0;
                return loolikelihood;
            };
            likeprofile_rms[ly][l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"rmslayer_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval").GetValue();
        }
    }

    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {   
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx; 
        auto map_filter = [lambda, l_idx, &h_sumrms_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hmean = h_sumrms_gauss[l_idx]->GetMean();
            auto hsigma = h_sumrms_gauss[l_idx]->GetRMS();
            new_val -= hmean>0 ? hmean : -hmean;
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto loglike = [l_idx, &h_sumrms_gauss_norm](double value) -> double 
        {
            double pdf = TMath::Gaus(value, h_sumrms_gauss_norm[l_idx]->GetMean(), h_sumrms_gauss_norm[l_idx]->GetRMS(), true);
            double loolikelihood = pdf>0 ? std::log(pdf) : 0;
            return loolikelihood;
        };
        likeprofile_sumrms[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"sumrms_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval").GetValue();
    }

    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {   
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx; 
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
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
            auto loglike = [l_idx, ly, &h_energyfrac_layer_gauss_norm](double value) -> double 
            {
                double pdf = TMath::Gaus(value, h_energyfrac_layer_gauss_norm[l_idx][ly]->GetMean(), h_energyfrac_layer_gauss_norm[l_idx][ly]->GetRMS(), true);
                double loolikelihood = pdf>0 ? std::log(pdf) : 0;
                return loolikelihood;
            };
            likeprofile_elf[ly][l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fraclayer_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval").GetValue();
        }
    }

    for (int l_idx=0; l_idx<=ell_lambda_values.num; ++l_idx)
    {   
        lambda = ell_lambda_values.start + ell_lambda_values.step*l_idx; 
        auto map_filter = [lambda, l_idx, &h_energyfrac_last_layer_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hmean = h_energyfrac_last_layer_gauss[l_idx]->GetMean();
            auto hsigma = h_energyfrac_last_layer_gauss[l_idx]->GetRMS();
            new_val -= hmean>0 ? hmean : -hmean;
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto loglike = [l_idx, &h_energyfrac_last_layer_gauss_norm](double value) -> double 
        {
            double pdf = TMath::Gaus(value, h_energyfrac_last_layer_gauss_norm[l_idx]->GetMean(), h_energyfrac_last_layer_gauss_norm[l_idx]->GetRMS(), true);
            double loolikelihood = pdf>0 ? std::log(pdf) : 0;
            return loolikelihood;
        };
        likeprofile_ell[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fraclastlayer_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval").GetValue();
    }
    
    for (int l_idx=0; l_idx<=xtrl_lambda_values.num; ++l_idx)
    {   
        lambda = xtrl_lambda_values.start + xtrl_lambda_values.step*l_idx; 
        auto map_filter = [lambda, l_idx, &h_xtrl_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hmean = h_xtrl_gauss[l_idx]->GetMean();
            auto hsigma = h_xtrl_gauss[l_idx]->GetRMS();
            new_val -= hmean>0 ? hmean : -hmean;
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto loglike = [l_idx, &h_xtrl_gauss_norm](double value) -> double 
        {
            double pdf = TMath::Gaus(value, h_xtrl_gauss_norm[l_idx]->GetMean(), h_xtrl_gauss_norm[l_idx]->GetRMS(), true);
            double loolikelihood = pdf>0 ? std::log(pdf) : 0;
            return loolikelihood;
        };
        likeprofile_xtrl[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"xtrl_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval").GetValue();
    }
    
    // Build lambda vectors
    std::vector<double> rms_lambdas (rms_lambda_values.num+1);
    std::vector<double> sumrms_lambdas (sumrms_lambda_values.num+1);
    std::vector<double> elf_lambdas (elf_lambda_values.num+1);
    std::vector<double> ell_lambdas (ell_lambda_values.num+1);
    std::vector<double> xtrl_lambdas (xtrl_lambda_values.num+1);

    for (int lidx=0; lidx<=rms_lambda_values.num; ++lidx)
        rms_lambdas[lidx] = rms_lambda_values.start + rms_lambda_values.step*lidx;
    for (int lidx=0; lidx<=sumrms_lambda_values.num; ++lidx)
        sumrms_lambdas[lidx] = sumrms_lambda_values.start + sumrms_lambda_values.step*lidx;
    for (int lidx=0; lidx<=elf_lambda_values.num; ++lidx)
        elf_lambdas[lidx] = elf_lambda_values.start + elf_lambda_values.step*lidx;
    for (int lidx=0; lidx<=ell_lambda_values.num; ++lidx)
        ell_lambdas[lidx] = ell_lambda_values.start + ell_lambda_values.step*lidx;
    for (int lidx=0; lidx<=xtrl_lambda_values.num; ++lidx)
        xtrl_lambdas[lidx] = xtrl_lambda_values.start + xtrl_lambda_values.step*lidx;
    
    // Build TGraphs
    std::vector<std::shared_ptr<TGraph>> rms_lambda_likelihood_profile (DAMPE_bgo_nLayers);
    std::vector<std::shared_ptr<TGraph>> elf_lambda_likelihood_profile (DAMPE_bgo_nLayers);
    
    for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
    {
        rms_lambda_likelihood_profile[ly] = std::make_shared<TGraph>(rms_lambda_values.num+1, &rms_lambdas[0], &likeprofile_rms[ly][0]);
        elf_lambda_likelihood_profile[ly] = std::make_shared<TGraph>(elf_lambda_values.num+1, &elf_lambdas[0], &likeprofile_elf[ly][0]);
        rms_lambda_likelihood_profile[ly]->SetName((std::string("gr_rms_layer_") + std::to_string(ly)).c_str());
        elf_lambda_likelihood_profile[ly]->SetName((std::string("gr_elf_layer_") + std::to_string(ly)).c_str());
    }

    std::shared_ptr<TGraph> sumrms_lambda_likelihood_profile = std::make_shared<TGraph>(sumrms_lambda_values.num+1, &sumrms_lambdas[0], &likeprofile_sumrms[0]);
    std::shared_ptr<TGraph> ell_lambda_likelihood_profile = std::make_shared<TGraph>(ell_lambda_values.num+1, &ell_lambdas[0], &likeprofile_ell[0]);
    std::shared_ptr<TGraph> xtrl_lambda_likelihood_profile = std::make_shared<TGraph>(xtrl_lambda_values.num+1, &xtrl_lambdas[0], &likeprofile_xtrl[0]);
    sumrms_lambda_likelihood_profile->SetName("gr_sumrms");
    ell_lambda_likelihood_profile->SetName("gr_ell");
    xtrl_lambda_likelihood_profile->SetName("gr_xtrl");

    output_file->mkdir((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/LProfiles")).c_str());
    output_file->cd((std::string("energybin_") + std::to_string(focus_energybin) + std::string("/LProfiles")).c_str());
    for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
    {
        rms_lambda_likelihood_profile[ly]->Write();
        elf_lambda_likelihood_profile[ly]->Write();
    }
    sumrms_lambda_likelihood_profile->Write();
    ell_lambda_likelihood_profile->Write();
    xtrl_lambda_likelihood_profile->Write();
    
    output_file->Close();   
    
    if (verbose)  std::cout << "\n\nOutput file has been written... [" << outputPath << "]\n";
}