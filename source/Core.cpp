#include "myHeader.h"

#include <vector>

void eCore(
            const std::string inputPath,
            const std::string outputPath,
            const bool verbose
            )
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

    readInputTree(inputPath,dataValues);

}