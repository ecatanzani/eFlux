#include "DAMPE_geo_structure.h"
#include "loglikelihood.h"
#include "histos.h"

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
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>>  h_rmsLayer_gauss (rms_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_sumrmsLayer_gauss (sumrms_lambda_values.num+1);
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>>  h_fracLayer_gauss (elf_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_fracLayer_ang_gauss (elf_ang_lambda_values.num+1);

    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_rmsLayer_gauss_norm (rms_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_sumrmsLayer_gauss_norm (sumrms_lambda_values.num+1);
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> h_fracLayer_gauss_norm (elf_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_fracLayer_ang_gauss_norm (elf_ang_lambda_values.num+1);
    
    double lambda;
    auto bin_filter = [focus_energybin](int energy_bin) -> bool { return energy_bin == (int)focus_energybin; };

    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {   
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx; 
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {   
            auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
            h_rmsLayer_gauss[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"rmsLayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {   
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx;   
        auto map_filter = [lambda](std::map<double, double> map_gauss) -> double { return map_gauss[lambda]; };
        h_sumrmsLayer_gauss[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"sumrmsLayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
    }

    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {   
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx; 
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
            h_fracLayer_gauss[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fracLayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=elf_ang_lambda_values.num; ++l_idx)
    {   
        lambda = elf_ang_lambda_values.start + elf_ang_lambda_values.step*l_idx; 
        auto map_filter = [lambda](std::map<double, double> map_gauss) -> double { return map_gauss[lambda]; };
        h_fracLayer_ang_gauss[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fracLayer_ang_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
    }
    
    output_file->mkdir("RMS");
    output_file->cd("RMS");
    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto h_name = std::string("h_rms_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);            
            h_rmsLayer_gauss[l_idx][ly]->SetName(h_name.c_str());
            h_rmsLayer_gauss[l_idx][ly]->GetXaxis()->SetTitle("RMS_{#lambda}");
            h_rmsLayer_gauss[l_idx][ly]->Write();
        }
    }

    output_file->mkdir("SumRMS");
    output_file->cd("SumRMS");
    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto h_name = std::string("h_sumrms_lambda_") + str_lambda;            
        h_sumrmsLayer_gauss[l_idx]->SetName(h_name.c_str());
        h_sumrmsLayer_gauss[l_idx]->GetXaxis()->SetTitle("SumRMS_{#lambda}");
        h_sumrmsLayer_gauss[l_idx]->Write();
    }
    
    output_file->mkdir("ELF");
    output_file->cd("ELF");
    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto h_name = std::string("h_fraclayer_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly); 
            h_fracLayer_gauss[l_idx][ly]->SetName(h_name.c_str());
            h_fracLayer_gauss[l_idx][ly]->GetXaxis()->SetTitle("ELF_{#lambda}");
            h_fracLayer_gauss[l_idx][ly]->Write();
        }
    }

    output_file->mkdir("ELFAng");
    output_file->cd("ELFAng");
    for (int l_idx=0; l_idx<=elf_ang_lambda_values.num; ++l_idx)
    {
        lambda = elf_ang_lambda_values.start + elf_ang_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto h_name = std::string("h_fraclayer_ang_lambda_") + str_lambda; 
        h_fracLayer_ang_gauss[l_idx]->SetName(h_name.c_str());
        h_fracLayer_ang_gauss[l_idx]->GetXaxis()->SetTitle("ELFang_{#lambda}");
        h_fracLayer_ang_gauss[l_idx]->Write();
    }

    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {   
            auto map_filter = [lambda, l_idx, ly, &h_rmsLayer_gauss](std::map<double, std::vector<double>> map_gauss) -> double 
            { 
                auto new_val = map_gauss[lambda][ly];
                auto hmean = h_rmsLayer_gauss[l_idx][ly]->GetMean();
                auto hsigma = h_rmsLayer_gauss[l_idx][ly]->GetRMS();
                new_val -= hmean>0 ? hmean : -hmean;
                if (hsigma) new_val /= hsigma;
                return new_val;
            };
            auto h_name = std::string("h_rms_norm_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
            auto h_title = std::string("Normalized RMS - #lambda ") + str_lambda + std::string(" - layer ") + std::to_string(ly);
            h_rmsLayer_gauss_norm[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                        .Define("mapval", map_filter, {"rmsLayer_gauss"})
                                                        .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto map_filter = [lambda, l_idx, &h_sumrmsLayer_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hmean = h_sumrmsLayer_gauss[l_idx]->GetMean();
            auto hsigma = h_sumrmsLayer_gauss[l_idx]->GetRMS();
            new_val -= hmean>0 ? hmean : -hmean;
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto h_name = std::string("h_sumrms_norm_lambda_") + str_lambda;
        auto h_title = std::string("Normalized SumRMS - #lambda ") + str_lambda;
        h_sumrmsLayer_gauss_norm[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("mapval", map_filter, {"sumrmsLayer_gauss"})
                                                    .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
    }

    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
        {   
            auto map_filter = [lambda, l_idx, ly, &h_fracLayer_gauss](std::map<double, std::vector<double>> map_gauss) -> double 
            { 
                auto new_val = map_gauss[lambda][ly];
                auto hmean = h_fracLayer_gauss[l_idx][ly]->GetMean();
                auto hsigma = h_fracLayer_gauss[l_idx][ly]->GetRMS();
                new_val -= hmean>0 ? hmean : -hmean;
                if (hsigma) new_val /= hsigma;
                return new_val;
            };
            auto h_name = std::string("h_fraclayer_norm_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);
            auto h_title = std::string("Normalized ELF - #lambda ") + str_lambda + std::string(" - layer ") + std::to_string(ly);
            h_fracLayer_gauss_norm[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                        .Define("mapval", map_filter, {"fracLayer_gauss"})
                                                        .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
        }
    }

    for (int l_idx=0; l_idx<=elf_ang_lambda_values.num; ++l_idx)
    {
        lambda = elf_ang_lambda_values.start + elf_ang_lambda_values.step*l_idx;
        auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
        auto map_filter = [lambda, l_idx, &h_fracLayer_ang_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hmean = h_fracLayer_ang_gauss[l_idx]->GetMean();
            auto hsigma = h_fracLayer_ang_gauss[l_idx]->GetRMS();
            new_val -= hmean>0 ? hmean : -hmean;
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto h_name = std::string("h_fraclayer_ang_norm_lambda_") + str_lambda;
        auto h_title = std::string("Normalized ELFang - #lambda ") + str_lambda;
        h_fracLayer_ang_gauss_norm[l_idx] = _data_fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("mapval", map_filter, {"fracLayer_ang_gauss"})
                                                    .Histo1D<double, double>({h_name.c_str(), h_title.c_str(), 100, -10, 10}, "mapval", "simu_energy_w_corr");
    }

    output_file->mkdir("RMS_norm");
    output_file->cd("RMS_norm");
    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            h_rmsLayer_gauss_norm[l_idx][ly]->Write();
            h_rmsLayer_gauss_norm[l_idx][ly]->GetXaxis()->SetTitle("RMS_{#lambda}");
        }
    
    output_file->mkdir("SumRMS_norm");
    output_file->cd("SumRMS_norm");
    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {
        h_sumrmsLayer_gauss_norm[l_idx]->Write();
        h_sumrmsLayer_gauss_norm[l_idx]->GetXaxis()->SetTitle("SumRMS_{#lambda}");
    }

    output_file->mkdir("ELF_norm");
    output_file->cd("ELF_norm");
    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            h_fracLayer_gauss_norm[l_idx][ly]->Write();
            h_fracLayer_gauss_norm[l_idx][ly]->GetXaxis()->SetTitle("ELF_{#lambda}");
        }

    output_file->mkdir("ELFang_norm");
    output_file->cd("ELFang_norm");
    for (int l_idx=0; l_idx<=elf_ang_lambda_values.num; ++l_idx)
    {
        h_fracLayer_ang_gauss_norm[l_idx]->Write();
        h_fracLayer_ang_gauss_norm[l_idx]->GetXaxis()->SetTitle("ELFang_{#lambda}");
    }

    if (verbose) std::cout << "\nBuilding RMS and ELF likelihood profiles...\n";
    std::vector<std::vector<ROOT::RDF::RResultPtr<double>>> likeprofile_rms (DAMPE_bgo_nLayers, std::vector<ROOT::RDF::RResultPtr<double>> (rms_lambda_values.num+1));
    std::vector<ROOT::RDF::RResultPtr<double>> likeprofile_sumrms (sumrms_lambda_values.num+1);
    std::vector<std::vector<ROOT::RDF::RResultPtr<double>>> likeprofile_elf (DAMPE_bgo_nLayers, std::vector<ROOT::RDF::RResultPtr<double>> (elf_lambda_values.num+1));
    std::vector<ROOT::RDF::RResultPtr<double>> likeprofile_elf_ang (elf_ang_lambda_values.num+1);
    
    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
    {   
        lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx; 
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto map_filter = [lambda, l_idx, ly, &h_rmsLayer_gauss](std::map<double, std::vector<double>> map_gauss) -> double 
            { 
                auto new_val = map_gauss[lambda][ly];
                auto hmean = h_rmsLayer_gauss[l_idx][ly]->GetMean();
                auto hsigma = h_rmsLayer_gauss[l_idx][ly]->GetRMS();
                new_val -= hmean>0 ? hmean : -hmean;
                if (hsigma) new_val /= hsigma;
                return new_val;
            };
            auto loglike = [l_idx, ly, &h_rmsLayer_gauss_norm](double value) -> double 
            {
                double pdf = TMath::Gaus(value, h_rmsLayer_gauss_norm[l_idx][ly]->GetMean(), h_rmsLayer_gauss_norm[l_idx][ly]->GetRMS(), true);
                double loolikelihood = pdf>0 ? std::log(pdf) : 0;
                return loolikelihood;
            };
            likeprofile_rms[ly][l_idx] = (_data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"rmsLayer_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval"));
        }
    }

    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
    {   
        lambda = sumrms_lambda_values.start + sumrms_lambda_values.step*l_idx; 
        auto map_filter = [lambda, l_idx, &h_sumrmsLayer_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hmean = h_sumrmsLayer_gauss[l_idx]->GetMean();
            auto hsigma = h_sumrmsLayer_gauss[l_idx]->GetRMS();
            new_val -= hmean>0 ? hmean : -hmean;
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto loglike = [l_idx, &h_sumrmsLayer_gauss_norm](double value) -> double 
        {
            double pdf = TMath::Gaus(value, h_sumrmsLayer_gauss_norm[l_idx]->GetMean(), h_sumrmsLayer_gauss_norm[l_idx]->GetRMS(), true);
            double loolikelihood = pdf>0 ? std::log(pdf) : 0;
            return loolikelihood;
        };
        likeprofile_sumrms[l_idx] = (_data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"sumrmsLayer_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval"));
    }

    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
    {   
        lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx; 
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            auto map_filter = [lambda, l_idx, ly, &h_fracLayer_gauss](std::map<double, std::vector<double>> map_gauss) -> double 
            { 
                auto new_val = map_gauss[lambda][ly];
                auto hmean = h_fracLayer_gauss[l_idx][ly]->GetMean();
                auto hsigma = h_fracLayer_gauss[l_idx][ly]->GetRMS();
                new_val -= hmean>0 ? hmean : -hmean;
                if (hsigma) new_val /= hsigma;
                return new_val;
            };
            auto loglike = [l_idx, ly, &h_fracLayer_gauss_norm](double value) -> double 
            {
                double pdf = TMath::Gaus(value, h_fracLayer_gauss_norm[l_idx][ly]->GetMean(), h_fracLayer_gauss_norm[l_idx][ly]->GetRMS(), true);
                double loolikelihood = pdf>0 ? std::log(pdf) : 0;
                return loolikelihood;
            };
            likeprofile_elf[ly][l_idx] = (_data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fracLayer_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval"));
        }
    }

    for (int l_idx=0; l_idx<=elf_ang_lambda_values.num; ++l_idx)
    {   
        lambda = elf_ang_lambda_values.start + elf_ang_lambda_values.step*l_idx; 
        auto map_filter = [lambda, l_idx, &h_fracLayer_ang_gauss](std::map<double, double> map_gauss) -> double 
        { 
            auto new_val = map_gauss[lambda];
            auto hmean = h_fracLayer_ang_gauss[l_idx]->GetMean();
            auto hsigma = h_fracLayer_ang_gauss[l_idx]->GetRMS();
            new_val -= hmean>0 ? hmean : -hmean;
            if (hsigma) new_val /= hsigma;
            return new_val;
        };
        auto loglike = [l_idx, &h_fracLayer_ang_gauss_norm](double value) -> double 
        {
            double pdf = TMath::Gaus(value, h_fracLayer_ang_gauss_norm[l_idx]->GetMean(), h_fracLayer_ang_gauss_norm[l_idx]->GetRMS(), true);
            double loolikelihood = pdf>0 ? std::log(pdf) : 0;
            return loolikelihood;
        };
        likeprofile_elf_ang[l_idx] = (_data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fracLayer_ang_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval"));
    }

    std::vector<std::vector<double>> flikeprofile_rms (DAMPE_bgo_nLayers, std::vector<double> (rms_lambda_values.num+1));
    std::vector<double> flikeprofile_sumrms (rms_lambda_values.num+1);
    std::vector<std::vector<double>> flikeprofile_elf (DAMPE_bgo_nLayers, std::vector<double> (elf_lambda_values.num+1));
    std::vector<double> flikeprofile_elf_ang (elf_lambda_values.num+1);

    for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            flikeprofile_rms[ly][l_idx] = *likeprofile_rms[ly][l_idx];
    
    for (int l_idx=0; l_idx<=sumrms_lambda_values.num; ++l_idx)
        flikeprofile_sumrms[l_idx] = *likeprofile_sumrms[l_idx];
        
    for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            flikeprofile_elf[ly][l_idx] = *likeprofile_elf[ly][l_idx];
    
    for (int l_idx=0; l_idx<=elf_ang_lambda_values.num; ++l_idx)
        flikeprofile_elf_ang[l_idx] = *likeprofile_elf_ang[l_idx];

    // Build lambda vectors
    std::vector<double> rms_lambdas (rms_lambda_values.num+1);
    std::vector<double> sumrms_lambdas (rms_lambda_values.num+1);
    std::vector<double> elf_lambdas (elf_lambda_values.num+1);
    std::vector<double> elf_lambdas_ang (elf_lambda_values.num+1);

    for (int lidx=0; lidx<=rms_lambda_values.num; ++lidx)
        rms_lambdas[lidx] = rms_lambda_values.start + rms_lambda_values.step*lidx;
    for (int lidx=0; lidx<=sumrms_lambda_values.num; ++lidx)
        sumrms_lambdas[lidx] = sumrms_lambda_values.start + sumrms_lambda_values.step*lidx;
    for (int lidx=0; lidx<=elf_lambda_values.num; ++lidx)
        elf_lambdas[lidx] = elf_lambda_values.start + elf_lambda_values.step*lidx;
    for (int lidx=0; lidx<=elf_ang_lambda_values.num; ++lidx)
        elf_lambdas_ang[lidx] = elf_ang_lambda_values.start + elf_ang_lambda_values.step*lidx;

    // Build TGraphs
    std::vector<std::shared_ptr<TGraph>> rms_lambda_likelihood_profile (DAMPE_bgo_nLayers);
    std::shared_ptr<TGraph> sumrms_lambda_likelihood_profile;
    std::vector<std::shared_ptr<TGraph>> elf_lambda_likelihood_profile (DAMPE_bgo_nLayers);
    std::shared_ptr<TGraph> elf_ang_lambda_likelihood_profile;

    for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
    {
        rms_lambda_likelihood_profile[ly] = std::make_shared<TGraph>(rms_lambda_values.num+1, &rms_lambdas[0], &flikeprofile_rms[ly][0]);
        elf_lambda_likelihood_profile[ly] = std::make_shared<TGraph>(elf_lambda_values.num+1, &elf_lambdas[0], &flikeprofile_elf[ly][0]);
        std::string rms_graph_name = std::string("gr_rms_layer_") + std::to_string(ly);
        std::string elf_graph_name = std::string("gr_elf_layer_") + std::to_string(ly);
        rms_lambda_likelihood_profile[ly]->SetName(rms_graph_name.c_str());
        elf_lambda_likelihood_profile[ly]->SetName(elf_graph_name.c_str());
    }

    sumrms_lambda_likelihood_profile = std::make_shared<TGraph>(sumrms_lambda_values.num+1, &sumrms_lambdas[0], &flikeprofile_sumrms[0]);
    elf_ang_lambda_likelihood_profile = std::make_shared<TGraph>(elf_ang_lambda_values.num+1, &elf_lambdas_ang[0], &flikeprofile_elf_ang[0]);
    sumrms_lambda_likelihood_profile->SetName("gr_sumrms");
    elf_ang_lambda_likelihood_profile->SetName("gr_elf_ang");

    output_file->mkdir("LProfiles/RMS");
    output_file->cd("LProfiles/RMS");
    for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        rms_lambda_likelihood_profile[ly]->Write();
    
    output_file->mkdir("LProfiles/SumRMS");
    output_file->cd("LProfiles/SumRMS");
    sumrms_lambda_likelihood_profile->Write();

    output_file->mkdir("LProfiles/ELF");
    output_file->cd("LProfiles/ELF");
    for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        elf_lambda_likelihood_profile[ly]->Write();
    
    output_file->mkdir("LProfiles/ELFang");
    output_file->cd("LProfiles/ELFang");
    elf_ang_lambda_likelihood_profile->Write();

    output_file->Close();   

    if (verbose)  std::cout << "\n\nOutput file has been written... [" << outputPath << "]\n";
}