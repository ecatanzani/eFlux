import os
import sys
import numpy as np
from ROOT import TTree, TFile, TBranch

def convertROOT(opts):
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
        

        data = TFile.Open(fNpath)
        for event in data.DmlNtup :
            energy.append(event.tt_bgoTotalEcorr_GeV/1000.)
            xtrl.append(event.tt_Xtrl)
    
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