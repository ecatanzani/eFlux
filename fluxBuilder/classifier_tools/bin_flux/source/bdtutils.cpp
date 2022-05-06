#include "bdtutils.h"

std::map<std::string, int> GetTMVAMethods(const std::string mymethod) {
    std::map<std::string, int> Use;

    Use["Cuts"] = 0;
    Use["CutsD"] = 0;
    Use["CutsPCA"] = 0;
    Use["CutsGA"] = 0;
    Use["CutsSA"] = 0;

    Use["Likelihood"] = 0;
    Use["LikelihoodD"] = 0;
    Use["LikelihoodPCA"] = 0;
    Use["LikelihoodKDE"] = 0;
    Use["LikelihoodMIX"] = 0;

    Use["PDERS"] = 0;
    Use["PDERSD"] = 0;
    Use["PDERSPCA"] = 0;
    Use["PDEFoam"] = 0;
    Use["PDEFoamBoost"] = 0;
    Use["KNN"] = 0;

    Use["LD"] = 0;
    Use["Fisher"] = 0;
    Use["FisherG"] = 0;
    Use["BoostedFisher"] = 0;
    Use["HMatrix"] = 0;

    Use["FDA_GA"] = 0;
    Use["FDA_SA"] = 0;
    Use["FDA_MC"] = 0;
    Use["FDA_MT"] = 0;
    Use["FDA_GAMT"] = 0;
    Use["FDA_MCMT"] = 0;

    Use["MLP"] = 0;
    Use["MLPBFGS"] = 0;
    Use["MLPBNN"] = 0;
    Use["CFMlpANN"] = 0;
    Use["TMlpANN"] = 0;

    Use["SVM"] = 0;

    Use["BDT"] = 0;
    Use["BDTG"] = 0;
    Use["BDTB"] = 0;
    Use["BDTD"] = 0;
    Use["BDTF"] = 0;

    Use["RuleFit"] = 0;

    auto linked_method = false;
    for (auto &&dmethod : Use)
        if (!strcmp(mymethod.c_str(), dmethod.first.c_str())) {
            dmethod.second = 1;
            linked_method = true;
            break;
        }

    if (!linked_method) {
        std::cerr << "\nERROR: No match found in TMVA default methods...\n\n";
        exit(100);
    }

    return Use;
}

void addVariableToReader(std::shared_ptr<TMVA::Reader> tmva_reader, bdt_vars &tmva_vars) {
    /*
    tmva_reader->AddVariable("rmslayer_norm_1", &tmva_vars.rmslayer_norm_1);
    tmva_reader->AddVariable("rmslayer_norm_2", &tmva_vars.rmslayer_norm_2);
    tmva_reader->AddVariable("rmslayer_norm_3", &tmva_vars.rmslayer_norm_3);
    tmva_reader->AddVariable("rmslayer_norm_4", &tmva_vars.rmslayer_norm_4);
    tmva_reader->AddVariable("rmslayer_norm_5", &tmva_vars.rmslayer_norm_5);
    tmva_reader->AddVariable("rmslayer_norm_6", &tmva_vars.rmslayer_norm_6);
    tmva_reader->AddVariable("rmslayer_norm_7", &tmva_vars.rmslayer_norm_7);
    tmva_reader->AddVariable("rmslayer_norm_8", &tmva_vars.rmslayer_norm_8);
    tmva_reader->AddVariable("rmslayer_norm_9", &tmva_vars.rmslayer_norm_9);
    tmva_reader->AddVariable("rmslayer_norm_10", &tmva_vars.rmslayer_norm_10);
    tmva_reader->AddVariable("rmslayer_norm_11", &tmva_vars.rmslayer_norm_11);
    tmva_reader->AddVariable("rmslayer_norm_12", &tmva_vars.rmslayer_norm_12);
    tmva_reader->AddVariable("rmslayer_norm_13", &tmva_vars.rmslayer_norm_13);
    tmva_reader->AddVariable("rmslayer_norm_14", &tmva_vars.rmslayer_norm_14);

    tmva_reader->AddVariable("fraclayer_norm_1", &tmva_vars.fraclayer_norm_1);
    tmva_reader->AddVariable("fraclayer_norm_2", &tmva_vars.fraclayer_norm_2);
    tmva_reader->AddVariable("fraclayer_norm_3", &tmva_vars.fraclayer_norm_3);
    tmva_reader->AddVariable("fraclayer_norm_4", &tmva_vars.fraclayer_norm_4);
    tmva_reader->AddVariable("fraclayer_norm_5", &tmva_vars.fraclayer_norm_5);
    tmva_reader->AddVariable("fraclayer_norm_6", &tmva_vars.fraclayer_norm_6);
    tmva_reader->AddVariable("fraclayer_norm_7", &tmva_vars.fraclayer_norm_7);
    tmva_reader->AddVariable("fraclayer_norm_8", &tmva_vars.fraclayer_norm_8);
    tmva_reader->AddVariable("fraclayer_norm_9", &tmva_vars.fraclayer_norm_9);
    tmva_reader->AddVariable("fraclayer_norm_10", &tmva_vars.fraclayer_norm_10);
    tmva_reader->AddVariable("fraclayer_norm_11", &tmva_vars.fraclayer_norm_11);
    tmva_reader->AddVariable("fraclayer_norm_12", &tmva_vars.fraclayer_norm_12);
    tmva_reader->AddVariable("fraclayer_norm_13", &tmva_vars.fraclayer_norm_13);
    tmva_reader->AddVariable("fraclayer_norm_14", &tmva_vars.fraclayer_norm_14);

    tmva_reader->AddVariable("sumrms_norm", &tmva_vars.sumrms_norm);
    tmva_reader->AddVariable("fraclastlayer_norm", &tmva_vars.fraclastlayer_norm);
    tmva_reader->AddVariable("xtrl_norm", &tmva_vars.xtrl_norm);
    tmva_reader->AddSpectator("xtrl", &tmva_vars.xtrl);
    */

    tmva_reader->AddVariable("rmslayer_1", &tmva_vars.rmslayer_norm_1);
    tmva_reader->AddVariable("rmslayer_2", &tmva_vars.rmslayer_norm_2);
    tmva_reader->AddVariable("rmslayer_3", &tmva_vars.rmslayer_norm_3);
    tmva_reader->AddVariable("rmslayer_4", &tmva_vars.rmslayer_norm_4);
    tmva_reader->AddVariable("rmslayer_5", &tmva_vars.rmslayer_norm_5);
    tmva_reader->AddVariable("rmslayer_6", &tmva_vars.rmslayer_norm_6);
    tmva_reader->AddVariable("rmslayer_7", &tmva_vars.rmslayer_norm_7);
    tmva_reader->AddVariable("rmslayer_8", &tmva_vars.rmslayer_norm_8);
    tmva_reader->AddVariable("rmslayer_9", &tmva_vars.rmslayer_norm_9);
    tmva_reader->AddVariable("rmslayer_10", &tmva_vars.rmslayer_norm_10);
    tmva_reader->AddVariable("rmslayer_11", &tmva_vars.rmslayer_norm_11);
    tmva_reader->AddVariable("rmslayer_12", &tmva_vars.rmslayer_norm_12);
    tmva_reader->AddVariable("rmslayer_13", &tmva_vars.rmslayer_norm_13);
    tmva_reader->AddVariable("rmslayer_14", &tmva_vars.rmslayer_norm_14);

    tmva_reader->AddVariable("fraclayer_1", &tmva_vars.fraclayer_norm_1);
    tmva_reader->AddVariable("fraclayer_2", &tmva_vars.fraclayer_norm_2);
    tmva_reader->AddVariable("fraclayer_3", &tmva_vars.fraclayer_norm_3);
    tmva_reader->AddVariable("fraclayer_4", &tmva_vars.fraclayer_norm_4);
    tmva_reader->AddVariable("fraclayer_5", &tmva_vars.fraclayer_norm_5);
    tmva_reader->AddVariable("fraclayer_6", &tmva_vars.fraclayer_norm_6);
    tmva_reader->AddVariable("fraclayer_7", &tmva_vars.fraclayer_norm_7);
    tmva_reader->AddVariable("fraclayer_8", &tmva_vars.fraclayer_norm_8);
    tmva_reader->AddVariable("fraclayer_9", &tmva_vars.fraclayer_norm_9);
    tmva_reader->AddVariable("fraclayer_10", &tmva_vars.fraclayer_norm_10);
    tmva_reader->AddVariable("fraclayer_11", &tmva_vars.fraclayer_norm_11);
    tmva_reader->AddVariable("fraclayer_12", &tmva_vars.fraclayer_norm_12);
    tmva_reader->AddVariable("fraclayer_13", &tmva_vars.fraclayer_norm_13);
    tmva_reader->AddVariable("fraclayer_14", &tmva_vars.fraclayer_norm_14);

    tmva_reader->AddVariable("sumRms", &tmva_vars.sumrms_norm);
    tmva_reader->AddVariable("fracLast", &tmva_vars.fraclastlayer_norm);
    tmva_reader->AddVariable("xtrl", &tmva_vars.xtrl_norm);
}

void sync_vars(const data_vars &vars, bdt_vars &tmva_vars) {
    tmva_vars.rmslayer_norm_1 = static_cast<float>(vars.rmslayer_norm_1);
    tmva_vars.rmslayer_norm_2 = static_cast<float>(vars.rmslayer_norm_2);
    tmva_vars.rmslayer_norm_3 = static_cast<float>(vars.rmslayer_norm_3);
    tmva_vars.rmslayer_norm_4 = static_cast<float>(vars.rmslayer_norm_4);
    tmva_vars.rmslayer_norm_5 = static_cast<float>(vars.rmslayer_norm_5);
    tmva_vars.rmslayer_norm_6 = static_cast<float>(vars.rmslayer_norm_6);
    tmva_vars.rmslayer_norm_7 = static_cast<float>(vars.rmslayer_norm_7);
    tmva_vars.rmslayer_norm_8 = static_cast<float>(vars.rmslayer_norm_8);
    tmva_vars.rmslayer_norm_9 = static_cast<float>(vars.rmslayer_norm_9);
    tmva_vars.rmslayer_norm_10 = static_cast<float>(vars.rmslayer_norm_10);
    tmva_vars.rmslayer_norm_11 = static_cast<float>(vars.rmslayer_norm_11);
    tmva_vars.rmslayer_norm_12 = static_cast<float>(vars.rmslayer_norm_12);
    tmva_vars.rmslayer_norm_13 = static_cast<float>(vars.rmslayer_norm_13);
    tmva_vars.rmslayer_norm_14 = static_cast<float>(vars.rmslayer_norm_14);

    tmva_vars.fraclayer_norm_1 = static_cast<float>(vars.fraclayer_norm_1);
    tmva_vars.fraclayer_norm_2 = static_cast<float>(vars.fraclayer_norm_2);
    tmva_vars.fraclayer_norm_3 = static_cast<float>(vars.fraclayer_norm_3);
    tmva_vars.fraclayer_norm_4 = static_cast<float>(vars.fraclayer_norm_4);
    tmva_vars.fraclayer_norm_5 = static_cast<float>(vars.fraclayer_norm_5);
    tmva_vars.fraclayer_norm_6 = static_cast<float>(vars.fraclayer_norm_6);
    tmva_vars.fraclayer_norm_7 = static_cast<float>(vars.fraclayer_norm_7);
    tmva_vars.fraclayer_norm_8 = static_cast<float>(vars.fraclayer_norm_8);
    tmva_vars.fraclayer_norm_9 = static_cast<float>(vars.fraclayer_norm_9);
    tmva_vars.fraclayer_norm_10 = static_cast<float>(vars.fraclayer_norm_10);
    tmva_vars.fraclayer_norm_11 = static_cast<float>(vars.fraclayer_norm_11);
    tmva_vars.fraclayer_norm_12 = static_cast<float>(vars.fraclayer_norm_12);
    tmva_vars.fraclayer_norm_13 = static_cast<float>(vars.fraclayer_norm_13);
    tmva_vars.fraclayer_norm_14 = static_cast<float>(vars.fraclayer_norm_14);

    tmva_vars.sumrms_norm = static_cast<float>(vars.sumrms_norm);
    tmva_vars.fraclastlayer_norm = static_cast<float>(vars.fraclastlayer_norm);
    tmva_vars.xtrl_norm = static_cast<float>(vars.xtrl_norm);
    tmva_vars.xtrl = static_cast<float>(vars.xtrl);
}

void linkTreeVariables(std::shared_ptr<TChain> evtch, data_vars &vars) {

    /*
    evtch->SetBranchAddress("rmslayer_norm_1", &vars.rmslayer_norm_1);
    evtch->SetBranchAddress("rmslayer_norm_2", &vars.rmslayer_norm_2);
    evtch->SetBranchAddress("rmslayer_norm_3", &vars.rmslayer_norm_3);
    evtch->SetBranchAddress("rmslayer_norm_4", &vars.rmslayer_norm_4);
    evtch->SetBranchAddress("rmslayer_norm_5", &vars.rmslayer_norm_5);
    evtch->SetBranchAddress("rmslayer_norm_6", &vars.rmslayer_norm_6);
    evtch->SetBranchAddress("rmslayer_norm_7", &vars.rmslayer_norm_7);
    evtch->SetBranchAddress("rmslayer_norm_8", &vars.rmslayer_norm_8);
    evtch->SetBranchAddress("rmslayer_norm_9", &vars.rmslayer_norm_9);
    evtch->SetBranchAddress("rmslayer_norm_10", &vars.rmslayer_norm_10);
    evtch->SetBranchAddress("rmslayer_norm_11", &vars.rmslayer_norm_11);
    evtch->SetBranchAddress("rmslayer_norm_12", &vars.rmslayer_norm_12);
    evtch->SetBranchAddress("rmslayer_norm_13", &vars.rmslayer_norm_13);
    evtch->SetBranchAddress("rmslayer_norm_14", &vars.rmslayer_norm_14);

    evtch->SetBranchAddress("fraclayer_norm_1", &vars.fraclayer_norm_1);
    evtch->SetBranchAddress("fraclayer_norm_2", &vars.fraclayer_norm_2);
    evtch->SetBranchAddress("fraclayer_norm_3", &vars.fraclayer_norm_3);
    evtch->SetBranchAddress("fraclayer_norm_4", &vars.fraclayer_norm_4);
    evtch->SetBranchAddress("fraclayer_norm_5", &vars.fraclayer_norm_5);
    evtch->SetBranchAddress("fraclayer_norm_6", &vars.fraclayer_norm_6);
    evtch->SetBranchAddress("fraclayer_norm_7", &vars.fraclayer_norm_7);
    evtch->SetBranchAddress("fraclayer_norm_8", &vars.fraclayer_norm_8);
    evtch->SetBranchAddress("fraclayer_norm_9", &vars.fraclayer_norm_9);
    evtch->SetBranchAddress("fraclayer_norm_10", &vars.fraclayer_norm_10);
    evtch->SetBranchAddress("fraclayer_norm_11", &vars.fraclayer_norm_11);
    evtch->SetBranchAddress("fraclayer_norm_12", &vars.fraclayer_norm_12);
    evtch->SetBranchAddress("fraclayer_norm_13", &vars.fraclayer_norm_13);
    evtch->SetBranchAddress("fraclayer_norm_14", &vars.fraclayer_norm_14);

    evtch->SetBranchAddress("sumrms_norm", &vars.sumrms_norm);
    evtch->SetBranchAddress("fraclastlayer_norm", &vars.fraclastlayer_norm);
    evtch->SetBranchAddress("xtrl_norm", &vars.xtrl_norm);
    evtch->SetBranchAddress("xtrl", &vars.xtrl);

    evtch->SetBranchAddress("energy_corr", &vars.evt_corr_energy);

    evtch->SetBranchAddress("energy_bin", &vars.energy_bin);
    */

    evtch->SetBranchAddress("rmslayer_1", &vars.rmslayer_norm_1);
    evtch->SetBranchAddress("rmslayer_2", &vars.rmslayer_norm_2);
    evtch->SetBranchAddress("rmslayer_3", &vars.rmslayer_norm_3);
    evtch->SetBranchAddress("rmslayer_4", &vars.rmslayer_norm_4);
    evtch->SetBranchAddress("rmslayer_5", &vars.rmslayer_norm_5);
    evtch->SetBranchAddress("rmslayer_6", &vars.rmslayer_norm_6);
    evtch->SetBranchAddress("rmslayer_7", &vars.rmslayer_norm_7);
    evtch->SetBranchAddress("rmslayer_8", &vars.rmslayer_norm_8);
    evtch->SetBranchAddress("rmslayer_9", &vars.rmslayer_norm_9);
    evtch->SetBranchAddress("rmslayer_10", &vars.rmslayer_norm_10);
    evtch->SetBranchAddress("rmslayer_11", &vars.rmslayer_norm_11);
    evtch->SetBranchAddress("rmslayer_12", &vars.rmslayer_norm_12);
    evtch->SetBranchAddress("rmslayer_13", &vars.rmslayer_norm_13);
    evtch->SetBranchAddress("rmslayer_14", &vars.rmslayer_norm_14);

    evtch->SetBranchAddress("fraclayer_1", &vars.fraclayer_norm_1);
    evtch->SetBranchAddress("fraclayer_2", &vars.fraclayer_norm_2);
    evtch->SetBranchAddress("fraclayer_3", &vars.fraclayer_norm_3);
    evtch->SetBranchAddress("fraclayer_4", &vars.fraclayer_norm_4);
    evtch->SetBranchAddress("fraclayer_5", &vars.fraclayer_norm_5);
    evtch->SetBranchAddress("fraclayer_6", &vars.fraclayer_norm_6);
    evtch->SetBranchAddress("fraclayer_7", &vars.fraclayer_norm_7);
    evtch->SetBranchAddress("fraclayer_8", &vars.fraclayer_norm_8);
    evtch->SetBranchAddress("fraclayer_9", &vars.fraclayer_norm_9);
    evtch->SetBranchAddress("fraclayer_10", &vars.fraclayer_norm_10);
    evtch->SetBranchAddress("fraclayer_11", &vars.fraclayer_norm_11);
    evtch->SetBranchAddress("fraclayer_12", &vars.fraclayer_norm_12);
    evtch->SetBranchAddress("fraclayer_13", &vars.fraclayer_norm_13);
    evtch->SetBranchAddress("fraclayer_14", &vars.fraclayer_norm_14);

    evtch->SetBranchAddress("sumRms", &vars.sumrms_norm);
    evtch->SetBranchAddress("fracLast", &vars.fraclastlayer_norm);
    evtch->SetBranchAddress("xtrl", &vars.xtrl_norm);

    evtch->SetBranchAddress("energy_corr", &vars.evt_corr_energy);

    evtch->SetBranchAddress("energy_bin", &vars.energy_bin);
}