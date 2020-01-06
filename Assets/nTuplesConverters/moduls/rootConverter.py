import os
import sys
from array import array
from ROOT import TChain,TTree,TFile,TBranch

def convertROOT(opts):
    
    # Import TTree branching module
    from handleTree import branchTree,fillTreeChain

    # Create output TFile
    outFile = TFile(opts.output,"RECREATE")
    if outFile.IsZombie():
        print('Error writing output file {}'.format(opts.output))
        sys.exit()

    # Create TChain to handle input ROOT files - "DmlNtup" is the TTree name in each ntuple ROOT file
    dChain = TChain("DmlNtup")

    # Add ROOT files to the chain
    for file in os.listdir(opts.input):
        fNpath = str(opts.input) + "/" + str(file)
        if opts.verbose:
            print('Adding {} to the chain...'.format(fNpath))
        dChain.AddFile(fNpath)
    
    nevents = dChain.GetEntries()
    if opts.verbose:
        print('Collected {} events ...'.format(nevents))
    
    # Create output TTree
    myTree = TTree("collectionTree","All Particle Tree")
    
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
    fillTreeChain(
                opts,
                dChain,
                nevents,
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
                kStep=1e+4
                )

    # Writing optput file
    outFile.Write()
    outFile.Close()
   