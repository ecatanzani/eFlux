import os
import sys
import numpy as np
from ROOT import TTree, TFile, TBranch

def convertPickle(opts):
    outFile = TFile(opts.output,"RECREATE")
    if outFile.IsZombie():
        print('Error writing output file {}'.format(opts.output))
        sys.exit()
    myTree = TTree("collectionTree","All Particle Tree")

    energy = []
    xtrl = []
    
    for file in os.listdir(opts.input):
        fNpath = str(opts.input) + "/" + str(file)
        if opts.verbose:
            print('Converting {}'.format(fNpath))
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

    myTree.Branch( 'eReco', t_energy, 'eReco/F' )
    myTree.Branch( 'xtrl', t_xtrl, 'xtrl/F' )

    myTree.Fill()

    outFile.Write()
    outFile.Close()
