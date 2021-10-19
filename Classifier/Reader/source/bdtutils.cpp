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