#!/usr/bin/env python

from ROOT import TTree, TFile, TBranch
import numpy as np
from array import array
import sys
import os

def main(ntuplesDir,outFilePath):
    outFile = TFile(outFilePath,"RECREATE")
    if outFile.IsZombie():
        print('Error writing output file {}'.format(outFilePath))
        sys.exit()
    myTree = TTree("collectionTree","All Particle Tree")

    energy = []
    xtrl = []
    
    for file in os.listdir(ntuplesDir):
        fNpath = str(ntuplesDir) + "/" + str(file)
        data = np.load(fNpath,'r')
        nRows = np.size(data,0)
        for evIdx in range(0,nRows):
            energy.append(data[evIdx,44]/1000.)
            xtrl.append(data[evIdx,49])
    
    t_energy = np.asarray(energy)
    t_xtrl = np.asarray(xtrl)
    nevents = len(energy)

    del energy
    del xtrl

    myTree.Branch( 'eReco', t_energy, 'eReco[nevents]/F' )
    myTree.Branch( 'xtrl', d, 'xtrl[nevents]/F' )

    myTree.Fill()

    outFile.Write()
    outFile.Close()

if __name__ == '__main__' :
    try:
        main(sys.argv[1],sys.argv[2])
    except AttributeError :
        print("--- ERROR IN : ", sys.argv[1],sys.argv[2])
        raise