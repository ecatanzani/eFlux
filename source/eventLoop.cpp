#include "myHeader.h"

void evLoop(TH2D &inputHisto,const std::string inputPath,const bool eClassifier)
{

    /* TTree variables

    double totalEnergy = 0;               vector position 0
    double totalEnergyCorr = 0;           vector position 1
    double xtrl = 0;                      vector position 2
    double satPositionX = 0;              vector position 3
    double satPositionY = 0;              vector position 4
    double satPositionZ = 0;              vector position 5
    double satVelocityX = 0;              vector position 6
    double satVelocityY = 0;              vector position 7
    double satVelocityZ = 0;              vector position 8

    */

    unsigned nData = 9;
    std::vector<double> dataValues(nData,0);

    //Create data TTree
    TTree* dTree = new TTree("collectionTree","Data Collection Tree");

    //Linking input TTree
    readInputTree(inputPath,dataValues,dTree);
   
    //Cleanup memory ...
    dTree->Delete();

}