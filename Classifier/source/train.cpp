#include "train.h"
#include "config.h"
#include "reader.h"
#include "train_utils.h"

#include "TChain.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/IMethod.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/DataLoader.h"

#include <memory>
#include <vector>

void Train(in_args input_args)
{
    // Read input data sets
    if (input_args.verbose)
        std::cout << "\n==> Reading signal training data [" << input_args.train_signal_input_list << "]" << std::endl;
    auto signal_train_tree = ReadTreeFromFile(input_args.train_signal_input_list, "trainSignal", input_args.verbose);
    if (input_args.verbose)
        std::cout << "\n==> Reading background training data [" << input_args.train_background_input_list << "]" << std::endl;
    auto background_train_tree = ReadTreeFromFile(input_args.train_background_input_list, "trainBackground", input_args.verbose);
    if (input_args.verbose)
        std::cout << "\n==> Reading signal test data [" << input_args.test_signal_input_list << "]" << std::endl;
    auto signal_test_tree = ReadTreeFromFile(input_args.test_signal_input_list, "testSignal", input_args.verbose);
    if (input_args.verbose)
        std::cout << "\n==> Reading background test data [" << input_args.test_background_input_list << "]" << std::endl;
    auto background_test_tree = ReadTreeFromFile(input_args.test_background_input_list, "testBackground", input_args.verbose);
    if (input_args.verbose)
        std::cout << "\n==> Writing TMVA output ROOT file [" << input_args.output_path << "]" << std::endl;

    if (input_args.verbose)
    {
        std::cout << "\n\n**** Data Statistics...\n";
        std::cout << "\nSignal TRAINING data set -> " << signal_train_tree->GetEntries() << " entries";
        std::cout << "\nSignal TEST data set -> " << signal_test_tree->GetEntries() << " entries";
        std::cout << "\nBackgound TRAINING data set -> " << background_train_tree->GetEntries() << " entries";
        std::cout << "\nBackgound TEST data set -> " << background_test_tree->GetEntries() << " entries";
        std::cout << "\n\n**************************\n";
    }

    // Loat TMVA library
    TMVA::Tools::Instance();

    // Parse config file
    std::unique_ptr<config> _config =
        std::make_unique<config>(
            input_args.config_dir,
            signal_train_tree,
            background_train_tree);

    if (input_args.verbose)
        _config->PrintVariableOptions();

    // Get TMVA methods
    auto methods_map = GetTMVAMethods(input_args.learning_method);

    // Create TMVA output file
    TFile *outfile = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (!outfile->IsOpen())
    {
        std::cout << "\n\nError writing output TTree: [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    if (input_args.verbose)
        std::cout << "\n==> Start TMVAClassification" << std::endl;

    std::shared_ptr<TMVA::Factory> factory = std::make_shared<TMVA::Factory>("TMVAClassification", outfile,
                                                                             "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
    std::shared_ptr<TMVA::DataLoader> dataloader = std::make_shared<TMVA::DataLoader>("dataset");

    // Define the input variables that shall be methods_mapd for the MVA training
    SetTMVAVariables(dataloader, _config->GetVariableOptions());

    // Set global event weights per tree
    double signalWeight = 1.0;
    double backgroundWeight = 1.0;

    // Add signal and backround Train/Test trees
    dataloader->AddSignalTree(signal_train_tree.get(), signalWeight, TMVA::Types::kTraining);
    dataloader->AddSignalTree(signal_test_tree.get(), signalWeight, TMVA::Types::kTesting);
    dataloader->AddBackgroundTree(background_train_tree.get(), backgroundWeight, TMVA::Types::kTraining);
    dataloader->AddBackgroundTree(background_test_tree.get(), backgroundWeight, TMVA::Types::kTesting);

    // Add events weights
    dataloader->SetSignalWeightExpression("simu_energy_w_corr");
    dataloader->SetBackgroundWeightExpression("simu_energy_w_corr");

    // Apply additional cuts on the signal and background samples
    TCut signal_cuts = "";
    TCut background_cuts = "";

#if 0
    SetTMVACuts(
        signal_cuts, 
        background_cuts, 
        input_args.verbose);
#endif

    auto loader_str =
        std::string("nTrain_Signal=") + std::to_string(_config->GetSignalTrainEvents()) +
        std::string(":nTrain_Background=") + std::to_string(_config->GetBackgroundTrainEvents()) +
        std::string(":nTest_Signal=") + std::to_string(_config->GetSignalTestEvents()) +
        std::string(":nTest_Background=") + std::to_string(_config->GetBackgroundTestEvents()) +
        std::string(":SplitMode=Random:NormMode=None:!V");
    dataloader->PrepareTrainingAndTestTree(signal_cuts, background_cuts, loader_str.c_str());

    // Book methods
    BookMethods(factory, dataloader, methods_map);

    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // ---- Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // ----- Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // Save the output
    outfile->Close();

    if (input_args.verbose)
        std::cout << "\n==> Classification done... Output file has been written [" << input_args.output_path << "]" << std::endl;
}