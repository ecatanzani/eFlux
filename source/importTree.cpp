#include "dataLoop.h"

std::shared_ptr<TTree> getDataTree(
    const std::string inputPath,
    std::vector<float> &dataValues)
{
    TFile inputTreeFile(inputPath.c_str(), "READ");
    if (!inputTreeFile.IsOpen())
    {
        std::cerr << "\n\nError opening input TTree: " << inputPath;
        exit(123);
    }
    auto myDataTree = std::shared_ptr<TTree>( static_cast <TTree *> (inputTreeFile.Get("collectionTree")) );
    inputTreeFile.Close();
    // Branch input Data Tree
    branchTree(*myDataTree, dataValues);

    return myDataTree;
}

void branchTree(TTree &myDataTree, std::vector<float> &dataValues)
{
    myDataTree.SetBranchAddress("eReco", &dataValues[0]);
    myDataTree.SetBranchAddress("eRecoCorr", &dataValues[1]);
    myDataTree.SetBranchAddress("xtrl", &dataValues[2]);
    myDataTree.SetBranchAddress("posizionX", &dataValues[3]);
    myDataTree.SetBranchAddress("posizionY", &dataValues[4]);
    myDataTree.SetBranchAddress("posizionZ", &dataValues[5]);
    myDataTree.SetBranchAddress("velocityX", &dataValues[6]);
    myDataTree.SetBranchAddress("velocityY", &dataValues[7]);
    myDataTree.SetBranchAddress("velocityZ", &dataValues[8]);
}