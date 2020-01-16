#!/usr/bin/env python

import sys
from datetime import datetime
from argparse import ArgumentParser

start=datetime.now()

def main(args=None):
    
    ### Parsing options
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="Random event extractor from TTree")

    parser.add_argument("-i","--input", type=str, dest='input', help='input TTree')
    parser.add_argument("-o","--output", type=str, dest='output', default="myRndmTree.root" , help='name of output root TTree')
    parser.add_argument("-v","--verbose", dest='verbose', default=False, action='store_true', help='run in high verbosity mode')
    parser.add_argument("-e","--events", type=int, dest='events', const=1000, nargs='?', help='Number of events to be extracted - default value 1000')

    opts = parser.parse_args(args)

    #Load analysis functions
    sys.path.append("moduls")
    from extractor import extractFromTree

    extractFromTree(opts)


if __name__ == '__main__' :
    main()
    print datetime.now()-start