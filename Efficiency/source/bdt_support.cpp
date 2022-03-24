#include "bdt_support.h"

#include <map>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iostream>

void check_bdt_learnign_method(const std::string method) {
    
    std::map<std::string, int> methods_map;

    methods_map["Cuts"] = 0;
    methods_map["CutsD"] = 0;
    methods_map["CutsPCA"] = 0;
    methods_map["CutsGA"] = 0;
    methods_map["CutsSA"] = 0;

    methods_map["Likelihood"] = 0;
    methods_map["LikelihoodD"] = 0;
    methods_map["LikelihoodPCA"] = 0;
    methods_map["LikelihoodKDE"] = 0;
    methods_map["LikelihoodMIX"] = 0;

    methods_map["PDERS"] = 0;
    methods_map["PDERSD"] = 0;
    methods_map["PDERSPCA"] = 0;
    methods_map["PDEFoam"] = 0;
    methods_map["PDEFoamBoost"] = 0;
    methods_map["KNN"] = 0;

    methods_map["LD"] = 0;
    methods_map["Fisher"] = 0;
    methods_map["FisherG"] = 0;
    methods_map["BoostedFisher"] = 0;
    methods_map["HMatrix"] = 0;

    methods_map["FDA_GA"] = 0;
    methods_map["FDA_SA"] = 0;
    methods_map["FDA_MC"] = 0;
    methods_map["FDA_MT"] = 0;
    methods_map["FDA_GAMT"] = 0;
    methods_map["FDA_MCMT"] = 0;

    methods_map["MLP"] = 0;
    methods_map["MLPBFGS"] = 0;
    methods_map["MLPBNN"] = 0;
    methods_map["CFMlpANN"] = 0;
    methods_map["TMlpANN"] = 0;

    methods_map["SVM"] = 0;

    methods_map["BDT"] = 0;
    methods_map["BDTG"] = 0;
    methods_map["BDTB"] = 0;
    methods_map["BDTD"] = 0;
    methods_map["BDTF"] = 0;

    methods_map["RuleFit"] = 0;

    auto linked_method = false;
    for (auto &&dmethod : methods_map)
        if (!strcmp(method.c_str(), dmethod.first.c_str())) {
            dmethod.second = 1;
            linked_method = true;
            break;
        }

    if (!linked_method) {
        std::cerr << "\nERROR: No match found in TMVA default methods...\n\n";
        exit(100);
    }
}

std::string parse_config_file(std::string bdt_config_file)
{
	std::ifstream input_file(bdt_config_file.c_str());
	if (!input_file.is_open())
	{
		std::cerr << "\nInput config file not found [" << bdt_config_file << "]\n\n";
		exit(100);
	}
	std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

const bdt_weights get_config_info(std::string parsed_config)
{
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);

    bdt_weights w;
	while (input_stream >> tmp_str)
	{
		if (!strcmp(tmp_str.c_str(), "low_energy_weights")) input_stream >> w.le_weights;
		if (!strcmp(tmp_str.c_str(), "medium_energy_weights")) input_stream >> w.me_weights;
		if (!strcmp(tmp_str.c_str(), "high_energy_weights")) input_stream >> w.he_weights;
	}

    return w;
}

void link_reader_vars(std::shared_ptr<TMVA::Reader> reader, classifier_vars &tmva_vars) {
    reader->AddVariable("rmslayer_norm_1", &tmva_vars.rmslayer_norm_1);
    reader->AddVariable("rmslayer_norm_2", &tmva_vars.rmslayer_norm_2);
    reader->AddVariable("rmslayer_norm_3", &tmva_vars.rmslayer_norm_3);
    reader->AddVariable("rmslayer_norm_4", &tmva_vars.rmslayer_norm_4);
    reader->AddVariable("rmslayer_norm_5", &tmva_vars.rmslayer_norm_5);
    reader->AddVariable("rmslayer_norm_6", &tmva_vars.rmslayer_norm_6);
    reader->AddVariable("rmslayer_norm_7", &tmva_vars.rmslayer_norm_7);
    reader->AddVariable("rmslayer_norm_8", &tmva_vars.rmslayer_norm_8);
    reader->AddVariable("rmslayer_norm_9", &tmva_vars.rmslayer_norm_9);
    reader->AddVariable("rmslayer_norm_10", &tmva_vars.rmslayer_norm_10);
    reader->AddVariable("rmslayer_norm_11", &tmva_vars.rmslayer_norm_11);
    reader->AddVariable("rmslayer_norm_12", &tmva_vars.rmslayer_norm_12);
    reader->AddVariable("rmslayer_norm_13", &tmva_vars.rmslayer_norm_13);
    reader->AddVariable("rmslayer_norm_14", &tmva_vars.rmslayer_norm_14);

    reader->AddVariable("fraclayer_norm_1", &tmva_vars.fraclayer_norm_1);
    reader->AddVariable("fraclayer_norm_2", &tmva_vars.fraclayer_norm_2);
    reader->AddVariable("fraclayer_norm_3", &tmva_vars.fraclayer_norm_3);
    reader->AddVariable("fraclayer_norm_4", &tmva_vars.fraclayer_norm_4);
    reader->AddVariable("fraclayer_norm_5", &tmva_vars.fraclayer_norm_5);
    reader->AddVariable("fraclayer_norm_6", &tmva_vars.fraclayer_norm_6);
    reader->AddVariable("fraclayer_norm_7", &tmva_vars.fraclayer_norm_7);
    reader->AddVariable("fraclayer_norm_8", &tmva_vars.fraclayer_norm_8);
    reader->AddVariable("fraclayer_norm_9", &tmva_vars.fraclayer_norm_9);
    reader->AddVariable("fraclayer_norm_10", &tmva_vars.fraclayer_norm_10);
    reader->AddVariable("fraclayer_norm_11", &tmva_vars.fraclayer_norm_11);
    reader->AddVariable("fraclayer_norm_12", &tmva_vars.fraclayer_norm_12);
    reader->AddVariable("fraclayer_norm_13", &tmva_vars.fraclayer_norm_13);
    reader->AddVariable("fraclayer_norm_14", &tmva_vars.fraclayer_norm_14);

    reader->AddVariable("sumrms_norm", &tmva_vars.sumrms_norm);
    reader->AddVariable("fraclastlayer_norm", &tmva_vars.fraclastlayer_norm);
    reader->AddVariable("xtrl_norm", &tmva_vars.xtrl_norm);
    reader->AddSpectator("xtrl", &tmva_vars.xtrl);
}

void bookMVA(
    std::shared_ptr<TMVA::Reader> reader, 
    const std::string weights, 
    const std::string method)
    {
        reader->BookMVA(method.c_str(), weights.c_str());
    }

void sync_vars(const bdt_vars &transformed_vars, classifier_vars &tmva_vars)
{
    tmva_vars.rmslayer_norm_1 = static_cast<float>(transformed_vars.rms[0]);
    tmva_vars.rmslayer_norm_2 = static_cast<float>(transformed_vars.rms[1]);
    tmva_vars.rmslayer_norm_3 = static_cast<float>(transformed_vars.rms[2]);
    tmva_vars.rmslayer_norm_4 = static_cast<float>(transformed_vars.rms[3]);
    tmva_vars.rmslayer_norm_5 = static_cast<float>(transformed_vars.rms[4]);
    tmva_vars.rmslayer_norm_6 = static_cast<float>(transformed_vars.rms[5]);
    tmva_vars.rmslayer_norm_7 = static_cast<float>(transformed_vars.rms[6]);
    tmva_vars.rmslayer_norm_8 = static_cast<float>(transformed_vars.rms[7]);
    tmva_vars.rmslayer_norm_9 = static_cast<float>(transformed_vars.rms[8]);
    tmva_vars.rmslayer_norm_10 = static_cast<float>(transformed_vars.rms[9]);
    tmva_vars.rmslayer_norm_11 = static_cast<float>(transformed_vars.rms[10]);
    tmva_vars.rmslayer_norm_12 = static_cast<float>(transformed_vars.rms[11]);
    tmva_vars.rmslayer_norm_13 = static_cast<float>(transformed_vars.rms[12]);
    tmva_vars.rmslayer_norm_14 = static_cast<float>(transformed_vars.rms[13]);

    tmva_vars.fraclayer_norm_1 = static_cast<float>(transformed_vars.fraclayer[0]);
    tmva_vars.fraclayer_norm_2 = static_cast<float>(transformed_vars.fraclayer[1]);
    tmva_vars.fraclayer_norm_3 = static_cast<float>(transformed_vars.fraclayer[2]);
    tmva_vars.fraclayer_norm_4 = static_cast<float>(transformed_vars.fraclayer[3]);
    tmva_vars.fraclayer_norm_5 = static_cast<float>(transformed_vars.fraclayer[4]);
    tmva_vars.fraclayer_norm_6 = static_cast<float>(transformed_vars.fraclayer[5]);
    tmva_vars.fraclayer_norm_7 = static_cast<float>(transformed_vars.fraclayer[6]);
    tmva_vars.fraclayer_norm_8 = static_cast<float>(transformed_vars.fraclayer[7]);
    tmva_vars.fraclayer_norm_9 = static_cast<float>(transformed_vars.fraclayer[8]);
    tmva_vars.fraclayer_norm_10 = static_cast<float>(transformed_vars.fraclayer[9]);
    tmva_vars.fraclayer_norm_11 = static_cast<float>(transformed_vars.fraclayer[10]);
    tmva_vars.fraclayer_norm_12 = static_cast<float>(transformed_vars.fraclayer[11]);
    tmva_vars.fraclayer_norm_13 = static_cast<float>(transformed_vars.fraclayer[12]);
    tmva_vars.fraclayer_norm_14 = static_cast<float>(transformed_vars.fraclayer[13]);

    tmva_vars.sumrms_norm = static_cast<float>(transformed_vars.sumrms);
    tmva_vars.fraclastlayer_norm = static_cast<float>(transformed_vars.fraclastlayer);
    tmva_vars.xtrl_norm = static_cast<float>(transformed_vars.xtrl);
    tmva_vars.xtrl = static_cast<float>(transformed_vars.xtrl_spectator);
}