#include "gaussianize.h"

#include <vector>
#include <cmath>
#include <map>

#include "DAMPE_geo_structure.h"

#include <ROOT/RDataFrame.hxx>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TFile.h"

void gaussianizeTMVAvars(
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> _config,
    std::shared_ptr<energy_config> _energy_config,
    std::shared_ptr<lambda_config> _lambda_config,
    const double _entries,
    const std::string outputPath,
    const bool _VERBOSE,
    const unsigned int threads,
    const bool _mc,
    const std::string tree_reg_path)
{
    // Enable multithreading
    ROOT::EnableImplicitMT(threads);
    // Create RDF
    ROOT::RDataFrame _data_fr(*evtch);
    // Extract the energy binning
    auto energy_binning = _energy_config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;
    // Create the RDFs
    double _gev = 0.001;
    auto GetEnergyBin = [=](double energy) -> int 
    { 
        int bin_idx=0;
        for (; bin_idx<energy_nbins-1; ++bin_idx)
            if (energy * _gev >= energy_binning[bin_idx] && energy * _gev < energy_binning[bin_idx+1])
                break;
        return bin_idx+1; 
    };
    // Create the RDF with energy bin and energy weight
    auto _fr_bin_patch = _data_fr.Define("energy_bin", GetEnergyBin, {"energy_corr"});
    if (_mc)
        _fr_bin_patch = _fr_bin_patch.Define("simu_energy_w_corr", [&energy_binning](const double simu_energy_w) -> double { return simu_energy_w * pow(energy_binning[0], 2); }, {"simu_energy_w"});
    else
        _fr_bin_patch = _fr_bin_patch.Define("simu_energy_w_corr", "1.");
    // Filtering events
    auto set_filter_min_energy = _energy_config->GetSetMinEvtEnergy();
    auto set_filter_max_energy = _energy_config->GetSetMaxEvtEnergy();
    auto _fr_preselected = _fr_bin_patch.Filter("evtfilter_all_cut==true")
                               .Filter([&set_filter_min_energy, &set_filter_max_energy, &_gev](const double energy) -> bool {
                                            auto status = false;
                                            if (energy*_gev >= set_filter_min_energy && energy*_gev <= set_filter_max_energy)
                                                status = true;
                                            else
                                                status = false;
                                            return status; }, {"energy_corr"});

    std::cout << "\n\n**** Filter statistics ****\n";
    std::cout << "***************************\n";
    std::cout << "\nPreselected events: " << *(_fr_preselected.Count());
    std::cout << "\n\n********************";

    // Load external TTree for TMVA variables regularization if needed
    auto load_gaus_reg_tree = [](std::string in_path) -> std::shared_ptr<TTreeReader>
    {
        TFile* regfile = TFile::Open(in_path.c_str(), "READ");
        if (regfile->IsZombie())
        {
            std::cerr << "\n\nError opening input gaus regularitazion TTree: [" << in_path << "]\n\n";
            exit(100);
        }
        std::shared_ptr<TTreeReader> _reader = std::make_shared<TTreeReader>("corrections_tree", regfile);
        return _reader;
    };
    std::shared_ptr<TTreeReader> _reader = !tree_reg_path.empty() ? load_gaus_reg_tree(tree_reg_path) : std::shared_ptr<TTreeReader>(nullptr);
    std::shared_ptr<TTreeReaderValue<std::vector<double>>> rms_lambda;
    std::shared_ptr<TTreeReaderValue<std::vector<double>>> rms_lambda_mean; 
    std::shared_ptr<TTreeReaderValue<std::vector<double>>> rms_lambda_sigma; 
    std::shared_ptr<TTreeReaderValue<std::vector<double>>> energyfraction_lambda;
    std::shared_ptr<TTreeReaderValue<std::vector<double>>> energyfraction_lambda_mean;
    std::shared_ptr<TTreeReaderValue<std::vector<double>>> energyfraction_lambda_sigma;
    if (_reader!=nullptr)
    {
        rms_lambda = std::make_shared<TTreeReaderValue<std::vector<double>>>(*_reader, "rms_lambda");
        rms_lambda_mean = std::make_shared<TTreeReaderValue<std::vector<double>>>(*_reader, "rms_lambda_mean");
        rms_lambda_sigma =  std::make_shared<TTreeReaderValue<std::vector<double>>>(*_reader, "rms_lambda_sigma");
        energyfraction_lambda = std::make_shared<TTreeReaderValue<std::vector<double>>>(*_reader, "energyfraction_lambda");
        energyfraction_lambda_mean = std::make_shared<TTreeReaderValue<std::vector<double>>>(*_reader, "energyfraction_lambda_mean");
        energyfraction_lambda_sigma = std::make_shared<TTreeReaderValue<std::vector<double>>>(*_reader, "energyfraction_lambda_sigma");
    }

    // Gaussianize TMVA variables
    auto rms_lambda_values = _lambda_config->GetRMSLambdaStruct();
    auto elf_lambda_values = _lambda_config->GetELFLambdaStruct();
    const bool iterative = true;
    const bool cycles = 100;
    if (_VERBOSE)
    {
        _lambda_config->PrintLambdaSettings();
        std::cout << "\nAnalysis running...\n";
    }

    auto gaussianize_rmslayer = [&rms_lambda_values, &iterative, &cycles](const std::vector<double> input_rmslayer) -> std::map<double, std::vector<double>>
    {
        auto gaussianize_elm = [](const std::vector<double> elm, const double lambda) -> std::vector<double> 
        {
            std::vector<double> elm_cp = elm;
            for (unsigned int idx=0; idx<elm.size(); ++idx)
                elm_cp[idx] = lambda ? (exp(lambda*elm[idx])-1)/lambda : elm[idx];
            return elm_cp;
        };
        
        std::map<double, std::vector<double>> rmsLayer_gauss;
        if (iterative)
        {
            for (unsigned int it_idx=0; it_idx<cycles; ++it_idx)
                for (int lambda_idx=0; lambda_idx<=rms_lambda_values.num; ++lambda_idx)
                {
                    double lambda = rms_lambda_values.start + lambda_idx*rms_lambda_values.step;
                    if (!it_idx) 
                        rmsLayer_gauss.insert(std::pair<double, std::vector<double>>(lambda, gaussianize_elm(input_rmslayer, lambda)));
                    else 
                        rmsLayer_gauss[lambda] = gaussianize_elm(rmsLayer_gauss[lambda], lambda);
                }
        }
        else
        {
            for (int lambda_idx=0; lambda_idx<=rms_lambda_values.num; ++lambda_idx)
            {
                double lambda = rms_lambda_values.start + lambda_idx*rms_lambda_values.step;
                rmsLayer_gauss.insert(std::pair<double, std::vector<double>>(lambda, gaussianize_elm(input_rmslayer, lambda)));
            }
        }

        return rmsLayer_gauss;
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

        std::map<double, std::vector<double>> fracLayer_gauss;
        if (iterative)
        {
            for (unsigned int it_idx=0; it_idx<cycles; ++it_idx)
                for (int lambda_idx=0; lambda_idx<=elf_lambda_values.num; ++lambda_idx)
                {
                    double lambda = elf_lambda_values.start + lambda_idx*elf_lambda_values.step;
                    if (!it_idx)
                        fracLayer_gauss.insert(std::pair<double, std::vector<double>>(lambda, gaussianize_elm(input_fraclayer, lambda)));
                    else
                        fracLayer_gauss[lambda] = gaussianize_elm(fracLayer_gauss[lambda], lambda);
                }
        }
        else
        {
            for (int lambda_idx=0; lambda_idx<=elf_lambda_values.num; ++lambda_idx)
            {
                double lambda = elf_lambda_values.start + lambda_idx*elf_lambda_values.step;
                fracLayer_gauss.insert(std::pair<double, std::vector<double>>(lambda, gaussianize_elm(input_fraclayer, lambda)));
            }
        }
        
        return fracLayer_gauss;
    };

    auto regularize_rmslayer = [&rms_lambda, &rms_lambda_mean, &rms_lambda_sigma, &_reader](const std::vector<double> rms, const int energy_idx) -> std::vector<double>
    {
        auto gaussianize_elm = [](
            const std::vector<double> elm, 
            const std::vector<double> lambda, 
            const std::vector<double> mean, 
            const::std::vector<double> sigma) -> std::vector<double> 
        {
            std::vector<double> elm_cp = elm;
            for (unsigned int idx=0; idx<elm.size(); ++idx)
            {
                elm_cp[idx] = lambda[idx] ? (exp(lambda[idx]*elm[idx])-1)/lambda[idx] : elm[idx];
                elm_cp[idx] -= mean[idx]>0 ? mean[idx] : -mean[idx];
                if (sigma[idx]) elm_cp[idx] /= sigma[idx];
            }
            return elm_cp;
        };
        _reader->SetEntry(energy_idx-1);
        return gaussianize_elm(rms, *(*rms_lambda), *(*rms_lambda_mean), *(*rms_lambda_sigma));
    };

    auto regularize_energyfractionlayer = [&energyfraction_lambda, &energyfraction_lambda_mean, &energyfraction_lambda_sigma, &_reader](const std::vector<double> energyfractionlayer, const int energy_idx) -> std::vector<double>
    {
        auto gaussianize_elm = [](
            const std::vector<double> elm, 
            const std::vector<double> lambda, 
            const std::vector<double> mean, 
            const::std::vector<double> sigma) -> std::vector<double> 
        {
            std::vector<double> elm_cp = elm;
            for (unsigned int idx=0; idx<elm.size(); ++idx)
            {
                elm_cp[idx] = lambda[idx] ? (exp(lambda[idx]*elm[idx])-1)/lambda[idx] : elm[idx];
                elm_cp[idx] -= mean[idx]>0 ? mean[idx] : -mean[idx];
                if (sigma[idx]) elm_cp[idx] /= sigma[idx];
            }
            return elm_cp;
        };
        _reader->SetEntry(energy_idx-1);
        return gaussianize_elm(energyfractionlayer, *(*energyfraction_lambda), *(*energyfraction_lambda_mean), *(*energyfraction_lambda_sigma));
    };

    auto _fr_preselected_gauss = tree_reg_path.empty() ? 
                                                _fr_preselected.Define("rmsLayer_gauss", gaussianize_rmslayer, {"rmsLayer"})
                                                .Define("fracLayer_gauss", gaussianize_fraclayer, {"fracLayer"}) : 
                                                _fr_preselected.Define("rmsLayer_gauss", regularize_rmslayer, {"rmsLayer", "energy_bin"})
                                                .Define("fracLayer_gauss", regularize_energyfractionlayer, {"fracLayer", "energy_bin"});

    _fr_preselected_gauss.Snapshot((std::string(evtch->GetName()) + std::string("_gauss")).c_str(), outputPath);
    
    if (_VERBOSE)
        std::cout << "\nOutput TFile has been written [" << outputPath << "]\n";
}