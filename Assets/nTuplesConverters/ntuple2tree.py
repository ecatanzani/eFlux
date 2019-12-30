#!/usr/bin/env python

import sys
from datetime import datetime
from argparse import ArgumentParser

start=datetime.now()

def main(args=None):
    
    ### Parsing options
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="nTuple to TTree converter")
    
    parser.add_argument("-i","--input",dest='input', help='ROOT nTuples directory')
    parser.add_argument("-o","--output",default="myTree.root", type=str, dest='outputFile', help='name of output root TTree')
    parser.add_argument("-v","--verbose", action='store_true', default=False, dest='verbose', help='run in high verbosity mode')
    parser.add_argument("-p","--pickle",action='store_true', default=False, dest='pickleFile', help='convert pickle files to TTree')
    parser.add_argument("-r","--root",action='store_true', default=False, dest='rootFile', help='convert root files to TTree')


    opts = parser.parse_args(args)

    #Load analysis functions
    sys.path.append("moduls")
    from pickleConverter import convertPickle
    from rootConverter import convertROOT

    if opts.pickle:
        if opts.verbose:
            print("Converting pickle numpy files to ROOT TTree...")
            convertPickle(opts)
    if opts.root:
        if opts.verbose:
            print("Converting ROOT files to ROOT TTree...")
            convertROOT(opts)


if __name__ == '__main__' :
    main()
    print datetime.now()-start
    