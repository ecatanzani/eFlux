#!/usr/bin/env python

import sys
from datetime import datetime
from argparse import ArgumentParser
from array import array
from tqdm import tqdm

from ROOT import TTree,TFile

start=datetime.now()

def fillTree(
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
                ):
    
    if opts.verbose:
        print("Filling output TTree...")
    
    #Reading input TFile
    inFile = TFile(opts.input,"READ")
    
    if not inFile.IsOpen():
        print('Error opening input TFile: {}'.format(opts.input))
        sys.exit()
    
    #Extracting input TTree
    dTree = inFile.Get("collectionTree")
    
    if opts.verbose:
        print('Collected {} events from input TTree'.format(dTree.GetEntries()))
        print('{} events will be random extracted'.format(opts.events))

    #Filling output TTree
    for idx,event in tqdm(enumerate(dTree)):
        t_totalEnergyCorr[0] = event.eRecoCorr
        t_totalEnergy[0] = event.eReco
        t_xtrl[0] = event.xtrl
        t_satPositionX[0] = event.posizionX
        t_satPositionY[0] = event.posizionY
        t_satPositionZ[0] = event.posizionZ
        t_satVelocityX[0] = event.velocityX
        t_satVelocityY[0] = event.velocityY
        t_satVelocityZ[0] = event.velocityZ
        myTree.Fill()
        if idx<(opts.events-1):
            idx +=1
        else:
            break
    
    #Close input TFile
    inFile.Close()

def branchTree(
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
                ):

    if opts.verbose:
        print("Branching output TTree...")
    myTree.Branch( 'eReco', t_totalEnergy, 'eReco/F' )
    myTree.Branch( 'eRecoCorr', t_totalEnergyCorr, 'eRecoCorr/F' )
    myTree.Branch( 'xtrl', t_xtrl, 'xtrl/F' )
    myTree.Branch( 'posizionX', t_satPositionX, 'posizionX/F' )
    myTree.Branch( 'posizionY', t_satPositionY, 'posizionY/F' )
    myTree.Branch( 'posizionZ', t_satPositionZ, 'posizionZ/F' )
    myTree.Branch( 'velocityX', t_satVelocityX, 'velocityX/F' )
    myTree.Branch( 'velocityY', t_satVelocityY, 'velocityY/F' )
    myTree.Branch( 'velocityZ', t_satVelocityZ, 'velocityZ/F' )

def extractFromTree(opts):

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
    fillTree(
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

    #Create output TFile
    outFile = TFile(opts.output,"RECREATE")
    if outFile.IsZombie():
        print('Error writing output file {}'.format(opts.output))
        sys.exit()
    
    #Add output TTree to the TFile
    myTree.Write()

    #Writing ouptput file
    outFile.Write()

    #Clone output TFile
    outFile.Close()


def main(args=None):
    
    ### Parsing options
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="Random event extractor from TTree")

    parser.add_argument("-i","--input", type=str, dest='input', help='input TTree')
    parser.add_argument("-o","--output", type=str, dest='output', default="myRndmTree.root" , help='name of output root TTree')
    parser.add_argument("-v","--verbose", dest='verbose', default=False, action='store_true', help='run in high verbosity mode')
    parser.add_argument("-e","--events", type=int, dest='events', const=1000, nargs='?', help='Number of events to be extracted - default value 1000')

    opts = parser.parse_args(args)

    extractFromTree(opts)


if __name__ == '__main__' :
    main()
    print datetime.now()-start