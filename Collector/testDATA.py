#!/usr/bin/env python3

from argparse import ArgumentParser
from ROOT import TFile

def TestROOTFile(path: str = "") -> bool:
    _tmp_file = TFile(path)
    if _tmp_file and not _tmp_file.IsOpen():
        return False
    elif _tmp_file and _tmp_file.IsOpen() and _tmp_file.IsZombie():
        _tmp_file.Close()
        return False
    elif _tmp_file and _tmp_file.IsOpen() and _tmp_file.TestBit(TFile.kRecovered):
        _tmp_file.Close()
        return False
    else:
        _tmp_file.Close()
        return True

def parseDataFiles(in_file: str = "", verbose: bool = False):
    zombie_files = []
    with open(in_file, "r") as _file:
        lines = _file.read().splitlines() 
    for line in lines:
        if verbose:
            print(f"Testing {line}")
        if not TestROOTFile(line):
            zombie_files.append(line)
            
    if zombie_files:
        print(f"\nHere the ZOMBIE data file list ({len(zombie_files)} elements / {len(lines)}):\n")
        with open("zombie_files.txt", "w") as zoutfile:
            for zfile in zombie_files:
                print(zfile)
                zoutfile.write(f"{zfile}\n")
            

def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="DATA test facility")

    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='Input DATA list')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')
    
    opts = parser.parse_args(args)

    if opts.input:
        parseDataFiles(opts.input, opts.verbose)
    else:
        print("Please select an input DATA list...")


if __name__ == "__main__":
    main()
