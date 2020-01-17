#include "myHeader.h"

void evLoop(
                    TFile &outFile,
                    TH1D &inputHisto,
                    const std::string inputPath,
                    const bool eClassifier,
                    const double xtrlCut
                )
{

    /* TTree variables

    double totalEnergy;               vector position 0
    double totalEnergyCorr;           vector position 1
    double xtrl;                      vector position 2
    double satPositionX;              vector position 3
    double satPositionY;              vector position 4
    double satPositionZ;              vector position 5
    double satVelocityX;              vector position 6
    double satVelocityY;              vector position 7
    double satVelocityZ;              vector position 8

    eClassifier True -> XTRL
    eClassifier False -> Manual Cuts

    */

    unsigned nData = 9;
    std::vector<double> dataValues(nData,0);

    //Create data TTree
    TTree* dTree = new TTree("collectionTree","Data Collection Tree");

    //Linking input TTree
    readInputTree(inputPath,dataValues,dTree);
    
    /*
        Telemetry information histos
            - Position of the satellite
            - Velocity of the satellite
    */

    TH1D satellitePosX("satellitePosX","satellite position X",1e+6,-8e+6,8e+6);
    TH1D satellitePosY("satellitePosY","satellite position Y",1e+6,-8e+6,8e+6);
    TH1D satellitePosZ("satellitePosZ","satellite position Z",1e+6,-8e+6,8e+6);

    TH1D satelliteVelX("satelliteVelX","satellite velocity X",1e+6,-8e+3,8e+3);
    TH1D satelliteVelY("satelliteVelY","satellite velocity Y",1e+6,-8e+3,8e+3);
    TH1D satelliteVelZ("satelliteVelZ","satellite velocity Z",1e+6,-8e+3,8e+3);


    for(unsigned int evIdx=0; evIdx<dTree->GetEntries(); ++evIdx)
    {
        dTree->GetEntry(evIdx);
        if(dataValues[2]<xtrlCut)
            inputHisto.Fill(dataValues[1]);
        
        satellitePosX.Fill(dataValues[3]);
        satellitePosY.Fill(dataValues[4]);
        satellitePosZ.Fill(dataValues[5]);

        satelliteVelX.Fill(dataValues[6]);
        satelliteVelY.Fill(dataValues[7]);
        satelliteVelZ.Fill(dataValues[8]);

    }

    //Cleanup memory ...
    dTree->Delete();

    //Write telemetry histos ...

    satellitePosX.Write();
    satellitePosY.Write();
    satellitePosZ.Write();

    satelliteVelX.Write();
    satelliteVelY.Write();
    satelliteVelZ.Write();

}