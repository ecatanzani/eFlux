#!/usr/bin/env python

import sys
from datetime import datetime
from argparse import ArgumentParser

start=datetime.now()

def main(args=None):
    
    ### Parsing options
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="nTuple to TTree converter")
    
    parser.add_argument("-i","--input", type=str, dest='input', help='ROOT nTuples directory')
    parser.add_argument("-o","--output", type=str, dest='output', default="myTree.root" , help='name of output root TTree')
    parser.add_argument("-v","--verbose", dest='verbose', default=False, action='store_true', help='run in high verbosity mode')
    parser.add_argument("-p","--pickle", dest='pickle', default=False, action='store_true', help='convert pickle files to TTree')
    parser.add_argument("-r","--root", dest='root', default=False, action='store_true', help='convert root files to TTree')
    #parser.add_argument("-d","--debug", dest='debug', action='store_const', const=1000, help='activate debug mode')
    parser.add_argument("-d","--debug", type=int, dest='debug', const=1000, nargs='?', help='activate debug mode')

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
    