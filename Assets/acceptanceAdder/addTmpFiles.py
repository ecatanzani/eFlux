import os
import sys
from datetime import datetime
from argparse import ArgumentParser
from ROOT import TFile, TH1D, TGraph
import math


start = datetime.now()


def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Acceptance Adder")

    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='Input condor jobs WD')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='Output ROOT TFile')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')
    parser.add_argument("-c", "--check", dest='check', default=False,
                        action='store_true', help='check tmp acceptance ROOT files')
    
    opts = parser.parse_args(args)

    # Load analysis functions
    sys.path.append("moduls")
    from adder import compute_final_histos
    from scanDirs import getListOfFiles

    good_dirs, skipped_dirs = getListOfFiles(opts.input)
    if opts.verbose:
        print('Found {} GOOD condor directories'.format(len(good_dirs)))

    if skipped_dirs:
        print('Found {} BAD condor directories...'.format(len(skipped_dirs)))
        for idx, elm in enumerate(skipped_dirs):
            print('Skipped {} directory: {}'.format(idx, elm))
    else:
        if opts.verbose:
            print('Found {} BAD condor directories...'.format(len(skipped_dirs)))

    if not opts.check:
        if opts.output:
            compute_final_histos(good_dirs, opts)
        else:
            print("!!! Missing output ROOT file... please specify a path")

if __name__ == "__main__":
    main()
    print(datetime.now()-start)
