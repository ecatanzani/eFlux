#include "myHeader.h"

/*
void readInputTree(const std::string inputPath,std::vector<float> &dataValues,TTree &dTree)
{
    TFile inputTree(inputPath.c_str(),"READ");
    if(!inputTree.IsOpen())
    {
        std::cerr << "\n\nError opening input TTree: " << inputPath;
        exit(123);
    }
    TTree *dataTree = (TTree*) inputTree.Get("collectionTree");
    new (&dTree) (TTree) (*(TTree*)dataTree->Clone("dataTree"));
    branchTree(dTree,dataValues);
}
*/

void branchTree(TTree &myDataTree,std::vector<float> &dataValues)
{
    myDataTree.SetBranchAddress("eReco",&dataValues[0]);
    myDataTree.SetBranchAddress("eRecoCorr",&dataValues[1]);
    myDataTree.SetBranchAddress("xtrl",&dataValues[2]);
    myDataTree.SetBranchAddress("posizionX",&dataValues[3]);
    myDataTree.SetBranchAddress("posizionY",&dataValues[4]);
    myDataTree.SetBranchAddress("posizionZ",&dataValues[5]);
    myDataTree.SetBranchAddress("velocityX",&dataValues[6]);
    myDataTree.SetBranchAddress("velocityY",&dataValues[7]);
    myDataTree.SetBranchAddress("velocityZ",&dataValues[8]);
}