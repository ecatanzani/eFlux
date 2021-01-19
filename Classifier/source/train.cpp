#include "train.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/IMethod.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/DataLoader.h"

int Train(in_args input_args)
{
    int status = 0;

    // Loat TMVA library
    TMVA::Tools::Instance();
    // Default MVA methods to be trained + tested
    std::map<std::string,int> Use;

    // --- Cut optimisation
    Use["Cuts"]            = 0;
    Use["CutsD"]           = 0;
    Use["CutsPCA"]         = 0;
    Use["CutsGA"]          = 0;
    Use["CutsSA"]          = 0;
    // 
    // --- 1-dimensional likelihood ("naive Bayes estimator")
    Use["Likelihood"]      = 0;
    Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
    Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
    Use["LikelihoodKDE"]   = 0;
    Use["LikelihoodMIX"]   = 0;
    //
    // --- Mutidimensional likelihood and Nearest-Neighbour methods
    Use["PDERS"]           = 0;
    Use["PDERSD"]          = 0;
    Use["PDERSPCA"]        = 0;
    Use["PDEFoam"]         = 0;
    Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
    Use["KNN"]             = 0; // k-nearest neighbour method
    //
    // --- Linear Discriminant Analysis
    Use["LD"]              = 0; // Linear Discriminant identical to Fisher
    Use["Fisher"]          = 0;
    Use["FisherG"]         = 0;
    Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
    Use["HMatrix"]         = 0;
    //
    // --- Function Discriminant analysis
    Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
    Use["FDA_SA"]          = 0;
    Use["FDA_MC"]          = 0;
    Use["FDA_MT"]          = 0;
    Use["FDA_GAMT"]        = 0;
    Use["FDA_MCMT"]        = 0;
    //
    // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
    Use["MLP"]             = 0; // Recommended ANN
    Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
    Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
    Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
    Use["TMlpANN"]         = 0; // ROOT's own ANN
    //
    // --- Support Vector Machine 
    Use["SVM"]             = 0;
    // 
    // --- Boosted Decision Trees
    Use["BDT"]             = 1; // uses Adaptive Boost
    Use["BDTG"]            = 1; // uses Gradient Boost
    Use["BDTB"]            = 0; // uses Bagging
    Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
    Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
    // 
    // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
    Use["RuleFit"]         = 0;
    // ---------------------------------------------------------------

    if (input_args.verbose)
        std::cout << "\n==> Reading signal training data set [" << input_args.train_signal_path << "]" << std::endl;
    auto signal_train_tree = ReadTreeFromFile(input_args.train_signal_path, "TreeSignalTrain");
    if (input_args.verbose)
        std::cout << "\n==> Reading background training data set [" << input_args.train_background_path << "]" << std::endl;
    auto background_train_tree = ReadTreeFromFile(input_args.train_background_path, "TreeBackgoundTrain");
    if (input_args.verbose)
        std::cout << "\n==> Reading signal test data set [" << input_args.test_signal_path << "]" << std::endl;
    auto signal_test_tree = ReadTreeFromFile(input_args.test_signal_path, "TreeSignalTest");
    if (input_args.verbose)
        std::cout << "\n==> Reading background test data set [" << input_args.test_background_path << "]" << std::endl;
    auto background_test_tree = ReadTreeFromFile(input_args.test_background_path, "TreeBackgoundTest");
    if (input_args.verbose)
        std::cout << "\n==> Writing TMVA output ROOT file [" << input_args.output_path << "]" << std::endl;
    TFile* outfile = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (!outfile->IsOpen())
    {
        std::cout << "\n\nError writing output TTree: [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    if (input_args.verbose)
        std::cout << "\n==> Start TMVAClassification" << std::endl;
    
    return status;
}

std::shared_ptr<TTree> ReadTreeFromFile(std::string tree_path, const char* tree_name)
{
    TFile *input = TFile::Open(tree_path.c_str(), "READ");
    if (!input->IsOpen())
    {
        std::cout << "\n\nError reading input TTree: [" << tree_path << "]\n\n";
        exit(100);
    }
    std::shared_ptr<TTree> tree = std::shared_ptr<TTree>(static_cast<TTree*>(input->Get(tree_name)));
    tree->SetDirectory(0);
    input->Close();
    return tree;
}