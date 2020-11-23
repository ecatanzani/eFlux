#!/usr/bin/env python3
import os
from argparse import ArgumentParser

def main(args=None):
    parser = ArgumentParser(
		usage="Usage: %(prog)s [options]", description="Create list of ROOT files from a folder")
    parser.add_argument("-i", "--input", type=str,
						dest='input', help='input folder')
    parser.add_argument("-o", "--output", type=str,
						dest='output', help='output file list')
    opts = parser.parse_args(args)
    
    _files = [ f for f in os.listdir(opts.input) if f.endswith('.root')]
    with open(opts.output, 'w') as _outfile:
        for _elm in _files:
            _outfile.write(f"{opts.input}/{_elm}\n")
    
if __name__ == "__main__":
    main()