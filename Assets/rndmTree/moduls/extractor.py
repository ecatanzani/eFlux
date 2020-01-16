from ROOT import TTree,TFile
from array import array

def extractFromTree(opts):

    from handleTree import branchTree,shuffleTree
    
    #Create output TFile
    outFile = TFile(opts.output,"RECREATE")
    if outFile.IsZombie():
        print('Error writing output file {}'.format(opts.output))
        sys.exit()

    # Create output TTree
    myTree = TTree("collectionTree","All Particle Tree - Rndm extracted")

    # Create final arrays for the TTree
    t_totalEnergy = array('f',[0.])
    t_totalEnergyCorr = array('f',[0.])
    t_xtrl = array('f',[0.])
    t_satPositionX = array('f',[0.])
    t_satPositionY = array('f',[0.])
    t_satPositionZ = array('f',[0.])
    t_satVelocityX = array('f',[0.])
    t_satVelocityY = array('f',[0.])
    t_satVelocityZ = array('f',[0.])

    # Branch TTree
    branchTree(
                opts,
                myTree,
                t_totalEnergyCorr,
                t_totalEnergy,
                t_xtrl,
                t_satPositionX,
                t_satPositionY,
                t_satPositionZ,
                t_satVelocityX,
                t_satVelocityY,
                t_satVelocityZ
                )

    # Filling the TTree
    shuffleTree(
                opts,
                myTree,
                t_totalEnergyCorr,
                t_totalEnergy,
                t_xtrl,
                t_satPositionX,
                t_satPositionY,
                t_satPositionZ,
                t_satVelocityX,
                t_satVelocityY,
                t_satVelocityZ,
                kStep=1e+2
                )

    #Writing ouptput file
    outFile.Write()

    #Clone output TFile
    outFile.Close()