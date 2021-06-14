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

#define single_bin true

void buildLogLikelihoodProfile(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    std::shared_ptr<lambda_config> _lambda_config,
    const double _entries,
    unsigned int focus_energybin,
    const std::string outputPath,
    const bool _VERBOSE,
    const unsigned int threads,
    const bool _mc)
{
    // Enable multithreading
    ROOT::EnableImplicitMT(threads);
    
    // Create RDF
    ROOT::RDataFrame _data_fr(*evtch);
    // Extract the energy binning
    auto energy_binning = _energy_config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;
    auto rms_lambda_values = _lambda_config->GetRMSLambdaStruct();
    auto elf_lambda_values = _lambda_config->GetELFLambdaStruct();
    if (_VERBOSE)
    {
        _lambda_config->PrintLambdaSettings();
        std::cout << "\nAnalysis running..." << std::endl;
    }

    if (single_bin)
    {
        TFile* output_file = TFile::Open(outputPath.c_str(), "RECREATE");
        if (output_file->IsZombie())
        {
            std::cerr << "\n\nError writing output file [" << outputPath << "]\n\n";
            exit(100);
        }

        if (_VERBOSE) std::cout << "\nBuilding RMS and LFS distributions..." << std::endl;
        std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>>  h_rmsLayer_gauss (rms_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
        std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>>  h_fracLayer_gauss (elf_lambda_values.num+1, std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers));
        
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

        for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
        {   
            lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx; 
            for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
                h_fracLayer_gauss[l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fracLayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
            }
        }
        
        output_file->mkdir("RMS");
        output_file->cd("RMS");
        for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
        {
            lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx;
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
                auto h_name = std::string("h_rms_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);            
                h_rmsLayer_gauss[l_idx][ly]->SetName(h_name.c_str());
                h_rmsLayer_gauss[l_idx][ly]->Write();
            }
        }
        
        output_file->mkdir("ELF");
        output_file->cd("ELF");
        for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
        {
            lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx;
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
                auto h_name = std::string("h_fraclayer_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly); 
                h_fracLayer_gauss[l_idx][ly]->SetName(h_name.c_str());
                h_fracLayer_gauss[l_idx][ly]->Write();
            }
        }

        if (_VERBOSE) std::cout << "\nBuilding RMS and ELF likelihood profiles...\n";
        std::vector<std::vector<ROOT::RDF::RResultPtr<double>>> likeprofile_rms (DAMPE_bgo_nLayers, std::vector<ROOT::RDF::RResultPtr<double>> (rms_lambda_values.num+1));
        std::vector<std::vector<ROOT::RDF::RResultPtr<double>>> likeprofile_elf (DAMPE_bgo_nLayers, std::vector<ROOT::RDF::RResultPtr<double>> (elf_lambda_values.num+1));
        
        for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
        {   
            lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx; 
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
                auto loglike = [l_idx, ly, &h_rmsLayer_gauss](double value) -> double 
                {
                    double pdf = TMath::Gaus(value, h_rmsLayer_gauss[l_idx][ly]->GetMean(), h_rmsLayer_gauss[l_idx][ly]->GetRMS(), true);
                    double loolikelihood = pdf>0 ? std::log(pdf) : 0;
                    return loolikelihood;
                };
                likeprofile_rms[ly][l_idx] = (_data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"rmsLayer_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval"));
            }
        }

        for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
        {   
            lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx; 
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
                auto loglike = [l_idx, ly, &h_fracLayer_gauss](double value) -> double 
                {
                    double pdf = TMath::Gaus(value, h_fracLayer_gauss[l_idx][ly]->GetMean(), h_fracLayer_gauss[l_idx][ly]->GetRMS(), true);
                    double loolikelihood = pdf>0 ? std::log(pdf) : 0;
                    return loolikelihood;
                };
                likeprofile_elf[ly][l_idx] = (_data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fracLayer_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval"));
            }
        }

        std::vector<std::vector<double>> flikeprofile_rms (DAMPE_bgo_nLayers, std::vector<double> (rms_lambda_values.num+1));
        std::vector<std::vector<double>> flikeprofile_elf (DAMPE_bgo_nLayers, std::vector<double> (elf_lambda_values.num+1));

        for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                flikeprofile_rms[ly][l_idx] = *likeprofile_rms[ly][l_idx];
            
        for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                flikeprofile_elf[ly][l_idx] = *likeprofile_elf[ly][l_idx];

        // Build lambda vectors
        std::vector<double> rms_lambdas (rms_lambda_values.num+1);
        std::vector<double> elf_lambdas (elf_lambda_values.num+1);

        for (int lidx=0; lidx<=rms_lambda_values.num; ++lidx)
            rms_lambdas[lidx] = rms_lambda_values.start + rms_lambda_values.step*lidx;
        for (int lidx=0; lidx<=elf_lambda_values.num; ++lidx)
            elf_lambdas[lidx] = elf_lambda_values.start + elf_lambda_values.step*lidx;

        // Build TGraphs
        std::vector<std::shared_ptr<TGraph>> rms_lambda_likelihood_profile (DAMPE_bgo_nLayers);
        std::vector<std::shared_ptr<TGraph>> elf_lambda_likelihood_profile (DAMPE_bgo_nLayers);

        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
        {
            rms_lambda_likelihood_profile[ly] = std::make_shared<TGraph>(rms_lambda_values.num+1, &rms_lambdas[0], &flikeprofile_rms[ly][0]);
            elf_lambda_likelihood_profile[ly] = std::make_shared<TGraph>(elf_lambda_values.num+1, &elf_lambdas[0], &flikeprofile_elf[ly][0]);
            std::string rms_graph_name = std::string("gr_rms_layer_") + std::to_string(ly);
            std::string elf_graph_name = std::string("gr_elf_layer_") + std::to_string(ly);
            rms_lambda_likelihood_profile[ly]->SetName(rms_graph_name.c_str());
            elf_lambda_likelihood_profile[ly]->SetName(elf_graph_name.c_str());
        }

        output_file->mkdir("LProfiles/RMS");
        output_file->cd("LProfiles/RMS");
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            rms_lambda_likelihood_profile[ly]->Write();
    
        output_file->mkdir("LProfiles/ELF");
        output_file->cd("LProfiles/ELF");
        for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            elf_lambda_likelihood_profile[ly]->Write();

        output_file->Close();
    }   
    else
    {
        TFile* output_file = TFile::Open(outputPath.c_str(), "RECREATE");
        if (output_file->IsZombie())
        {
            std::cerr << "\n\nError writing output file [" << outputPath << "]\n\n";
            exit(100);
        }

        if (_VERBOSE) std::cout << "\nBuilding RMS and LFS distributions..." << std::endl;
        std::vector<std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>>>  h_rmsLayer_gauss (energy_nbins, std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> (rms_lambda_values.num+1));
        std::vector<std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>>>  h_fracLayer_gauss (energy_nbins, std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> (elf_lambda_values.num+1));

        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
        {
            for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
                h_rmsLayer_gauss[bin_idx-1][l_idx] = std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers);
            for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
                h_fracLayer_gauss[bin_idx-1][l_idx] = std::vector<ROOT::RDF::RResultPtr<TH1D>> (DAMPE_bgo_nLayers);
        }

        double lambda;
        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
        {
            auto bin_filter = [bin_idx](int energy_bin) -> bool { return energy_bin == bin_idx; };
            for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
            {   
                lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx; 
                for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
                {   
                    auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
                    h_rmsLayer_gauss[bin_idx-1][l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"rmsLayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
                }
            }

            for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
            {   
                lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx; 
                for (int ly=0; ly < DAMPE_bgo_nLayers; ++ly)
                {
                    auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
                    h_fracLayer_gauss[bin_idx-1][l_idx][ly] = _data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fracLayer_gauss"}).Histo1D<double, double>("mapval", "simu_energy_w_corr");
                }
            }
        }
        
        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
        {
            output_file->mkdir((std::string("RMS/energybin_") + std::to_string(bin_idx)).c_str());
            output_file->cd((std::string("RMS/energybin_") + std::to_string(bin_idx)).c_str());
            for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
            {
                lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx;
                for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                {
                    auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
                    auto h_name = std::string("h_rms_energybin_") + std::to_string(bin_idx) + std::string("_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly);            
                    h_rmsLayer_gauss[bin_idx-1][l_idx][ly]->SetName(h_name.c_str());
                    h_rmsLayer_gauss[bin_idx-1][l_idx][ly]->Write();
                }
            }
            
            output_file->mkdir((std::string("ELF/energybin_") + std::to_string(bin_idx)).c_str());
            output_file->cd((std::string("ELF/energybin_") + std::to_string(bin_idx)).c_str());
            for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
            {
                lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx;
                for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                {
                    auto str_lambda = lambda<0 ? std::string("neg_") + std::to_string(std::abs(lambda)) : std::to_string(lambda);
                    auto h_name = std::string("h_fraclayer_energybin_") + std::to_string(bin_idx) + std::string("_lambda_") + str_lambda + std::string("_layer_") + std::to_string(ly); 
                    h_fracLayer_gauss[bin_idx-1][l_idx][ly]->SetName(h_name.c_str());
                    h_fracLayer_gauss[bin_idx-1][l_idx][ly]->Write();
                }
            }
        }
        
    #if 0
        std::vector<std::vector<std::vector<std::shared_ptr<TF1>>>> rms_fit_functions (energy_nbins, std::vector<std::vector<std::shared_ptr<TF1>>> (DAMPE_bgo_nLayers));
        std::vector<std::vector<std::vector<std::shared_ptr<TF1>>>> elf_fit_functions (energy_nbins, std::vector<std::vector<std::shared_ptr<TF1>>> (DAMPE_bgo_nLayers));
        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                rms_fit_functions[bin_idx-1][ly] = std::vector<std::shared_ptr<TF1>> (rms_lambda_values.num+1);
                elf_fit_functions[bin_idx-1][ly] = std::vector<std::shared_ptr<TF1>> (elf_lambda_values.num+1);
            }

        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
        {
            for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
                for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                {
                    auto xmin = h_rmsLayer_gauss[bin_idx-1][l_idx][ly]->GetXaxis()->GetXmin();
                    auto xmax = h_rmsLayer_gauss[bin_idx-1][l_idx][ly]->GetXaxis()->GetXmax();
                    std::string func_name = std::string("fitfunc_rms_energybin_") + std::to_string(bin_idx) + std::string("_ly_") + std::to_string(ly) + std::string("_lidx_") + std::to_string(l_idx);
                    rms_fit_functions[bin_idx-1][ly][l_idx] = std::make_shared<TF1>(func_name.c_str(), "gaus", xmin, xmax);
                    rms_fit_functions[bin_idx-1][ly][l_idx]->SetParameters(1, h_rmsLayer_gauss[bin_idx-1][l_idx][ly]->GetMean(), h_rmsLayer_gauss[bin_idx-1][l_idx][ly]->GetRMS());
                    h_rmsLayer_gauss[bin_idx-1][l_idx][ly]->Fit(func_name.c_str(), "QW");   
                }
            
            for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
                for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                {
                    auto xmin = h_fracLayer_gauss[bin_idx-1][l_idx][ly]->GetXaxis()->GetXmin();
                    auto xmax = h_fracLayer_gauss[bin_idx-1][l_idx][ly]->GetXaxis()->GetXmax();
                    std::string func_name = std::string("fitfunc_elf_energybin_") + std::to_string(bin_idx) + std::string("_ly_") + std::to_string(ly) + std::string("_lidx_") + std::to_string(l_idx);
                    elf_fit_functions[bin_idx-1][ly][l_idx] = std::make_shared<TF1>(func_name.c_str(), "gaus", xmin, xmax);
                    elf_fit_functions[bin_idx-1][ly][l_idx]->SetParameters(1, h_fracLayer_gauss[bin_idx-1][l_idx][ly]->GetMean(), h_fracLayer_gauss[bin_idx-1][l_idx][ly]->GetRMS());
                    h_fracLayer_gauss[bin_idx-1][l_idx][ly]->Fit(func_name.c_str(), "QW");
                }
        }
    #endif

        if (_VERBOSE) std::cout << "\nBuilding RMS and ELF likelihood profiles...\n";
        std::vector<std::vector<std::vector<ROOT::RDF::RResultPtr<double>>>> likeprofile_rms (energy_nbins, std::vector<std::vector<ROOT::RDF::RResultPtr<double>>> (DAMPE_bgo_nLayers));
        std::vector<std::vector<std::vector<ROOT::RDF::RResultPtr<double>>>> likeprofile_elf (energy_nbins, std::vector<std::vector<ROOT::RDF::RResultPtr<double>>> (DAMPE_bgo_nLayers));

        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                likeprofile_rms[bin_idx-1][ly] = std::vector<ROOT::RDF::RResultPtr<double>> (rms_lambda_values.num+1);
                likeprofile_elf[bin_idx-1][ly] = std::vector<ROOT::RDF::RResultPtr<double>> (elf_lambda_values.num+1);
            }
        
        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
        {
            auto bin_filter = [bin_idx](int energy_bin) -> bool { return energy_bin == bin_idx; };
            for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
            {   
                lambda = rms_lambda_values.start + rms_lambda_values.step*l_idx; 
                for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                {
                    auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
                    auto loglike = [bin_idx, l_idx, ly, &h_rmsLayer_gauss](double value) -> double 
                    {
                        double pdf = TMath::Gaus(value, h_rmsLayer_gauss[bin_idx-1][l_idx][ly]->GetMean(), h_rmsLayer_gauss[bin_idx-1][l_idx][ly]->GetRMS(), true);
                        double loolikelihood = pdf>0 ? std::log(pdf) : 0;
                        return loolikelihood;
                    };
                    likeprofile_rms[bin_idx-1][ly][l_idx] = (_data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"rmsLayer_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval"));
                }
            }

            for (int l_idx=0; l_idx<=elf_lambda_values.num; ++l_idx)
            {   
                lambda = elf_lambda_values.start + elf_lambda_values.step*l_idx; 
                for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                {
                    auto map_filter = [lambda, ly](std::map<double, std::vector<double>> map_gauss) -> double { return map_gauss[lambda][ly]; };
                    auto loglike = [bin_idx, l_idx, ly, &h_fracLayer_gauss](double value) -> double 
                    {
                        double pdf = TMath::Gaus(value, h_fracLayer_gauss[bin_idx-1][l_idx][ly]->GetMean(), h_fracLayer_gauss[bin_idx-1][l_idx][ly]->GetRMS(), true);
                        double loolikelihood = pdf>0 ? std::log(pdf) : 0;
                        return loolikelihood;
                    };
                    likeprofile_elf[bin_idx-1][ly][l_idx] = (_data_fr.Filter(bin_filter, {"energy_bin"}).Define("mapval", map_filter, {"fracLayer_gauss"}).Define("likeval", loglike, {"mapval"}).Sum<double>("likeval"));
                }
            }
        }
        
        std::vector<std::vector<std::vector<double>>> flikeprofile_rms (energy_nbins, std::vector<std::vector<double>> (DAMPE_bgo_nLayers));
        std::vector<std::vector<std::vector<double>>> flikeprofile_elf (energy_nbins, std::vector<std::vector<double>> (DAMPE_bgo_nLayers));
        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
        {
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                flikeprofile_rms[bin_idx-1][ly] = std::vector<double> (rms_lambda_values.num+1);
                flikeprofile_elf[bin_idx-1][ly] = std::vector<double> (elf_lambda_values.num+1);
            }

            for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
                for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                    flikeprofile_rms[bin_idx-1][ly][l_idx] = *likeprofile_rms[bin_idx-1][ly][l_idx];
            
            for (int l_idx=0; l_idx<=rms_lambda_values.num; ++l_idx)
                for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                    flikeprofile_elf[bin_idx-1][ly][l_idx] = *likeprofile_elf[bin_idx-1][ly][l_idx];
                
        }
        
        // Build lambda vectors
        std::vector<double> rms_lambdas (rms_lambda_values.num+1);
        std::vector<double> elf_lambdas (elf_lambda_values.num+1);

        for (int lidx=0; lidx<=rms_lambda_values.num; ++lidx)
            rms_lambdas[lidx] = rms_lambda_values.start + rms_lambda_values.step*lidx;
        for (int lidx=0; lidx<=elf_lambda_values.num; ++lidx)
            elf_lambdas[lidx] = elf_lambda_values.start + elf_lambda_values.step*lidx;
        
        // Build TGraphs
        std::vector<std::vector<std::shared_ptr<TGraph>>> rms_lambda_likelihood_profile (energy_nbins, std::vector<std::shared_ptr<TGraph>> (DAMPE_bgo_nLayers));
        std::vector<std::vector<std::shared_ptr<TGraph>>> elf_lambda_likelihood_profile (energy_nbins, std::vector<std::shared_ptr<TGraph>> (DAMPE_bgo_nLayers));

        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
        {
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
            {
                rms_lambda_likelihood_profile[bin_idx-1][ly] = std::make_shared<TGraph>(rms_lambda_values.num+1, &rms_lambdas[0], &flikeprofile_rms[bin_idx-1][ly][0]);
                elf_lambda_likelihood_profile[bin_idx-1][ly] = std::make_shared<TGraph>(elf_lambda_values.num+1, &elf_lambdas[0], &flikeprofile_elf[bin_idx-1][ly][0]);
                std::string rms_graph_name = std::string("gr_rms_energybin_") + std::to_string(bin_idx) + std::string("_layer_") + std::to_string(ly);
                std::string elf_graph_name = std::string("gr_elf_energybin_") + std::to_string(bin_idx) + std::string("_layer_") + std::to_string(ly);
                rms_lambda_likelihood_profile[bin_idx-1][ly]->SetName(rms_graph_name.c_str());
                elf_lambda_likelihood_profile[bin_idx-1][ly]->SetName(elf_graph_name.c_str());
            }
        }

        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx)
        {
            output_file->mkdir((std::string("LProfiles/RMS/energybin_") + std::to_string(bin_idx)).c_str());
            output_file->cd((std::string("LProfiles/RMS/energybin_") + std::to_string(bin_idx)).c_str());
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                rms_lambda_likelihood_profile[bin_idx-1][ly]->Write();
        
            output_file->mkdir((std::string("LProfiles/ELF/energybin_") + std::to_string(bin_idx)).c_str());
            output_file->cd((std::string("LProfiles/ELF/energybin_") + std::to_string(bin_idx)).c_str());
            for (int ly = 0; ly < DAMPE_bgo_nLayers; ++ly)
                elf_lambda_likelihood_profile[bin_idx-1][ly]->Write();
        }

        output_file->Close();
    }

    if (_VERBOSE)  std::cout << "\n\nOutput file has been written... [" << outputPath << "]\n";
}