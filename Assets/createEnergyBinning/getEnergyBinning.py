#!/usr/bin/env python3
from datetime import datetime
from argparse import ArgumentParser
import numpy as np
import sys
import math
import json

start = datetime.now()


def createBinning(opts, pars):
    best_binning = []
    min_power = math.log10(float(pars["eMin"]))/math.log10(10)
    max_power = math.log10(float(pars["eMax"]))/math.log10(10)

    if not pars["junctions"]:
        best_binning.append(np.logspace(min_power, max_power, pars["bins"]+1))
    else:
        max_error = 10
        min_bins = pars["bins"]-5
        max_bins = pars["bins"]+5
        if min_bins < 0 | max_bins<0:
            print("Error defining number of max or min bins... exit now")
            sys.exit()
        for error in range(max_error, 0, -1):
            best_error = []
            for idx in range(0, len(pars["junctions"])):
                best_error.append(error)
            for bins in range(min_bins, max_bins):
                binning = np.logspace(min_power, max_power, bins+1)
                found = [ False for _ in range(len(pars["junctions"]))]
                for bin_value in binning:
                    for jIdx in range(len(pars["junctions"])):
                        tmp_diff = np.absolute(bin_value-pars["junctions"][jIdx])
                        if (tmp_diff) <= best_error[jIdx]:
                            found[jIdx] = True
                            if (tmp_diff) < best_error[jIdx]:
                                best_error[jIdx] = np.absolute(tmp_diff)
                allJunctions = True
                for elm in found:
                    allJunctions *= elm
                if allJunctions:
                    best_binning.append(binning)

    if opts.verbose:
        if len(best_binning):
            print("Binning found!")
            print("Number of bins: ", len(best_binning[-1])-1)
            print("****\n")
            print('{}\n'.format(best_binning[-1]))

    if len(best_binning):
        return best_binning[-1]
    else:
        return best_binning


def main(args=None):

    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Binning Finder")

    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='json data storage')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='Output binning text file')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)

    # Parsing config file
    configFile = "config.json"
    if opts.input:
        configFile = opts.input

    with open(configFile, 'rb') as set_:
            pars = json.load(set_)
    
    if opts.verbose:
        print("\n**** Binning settings:")
        print('eMin: {}'.format(pars["eMin"]))
        print('eMax: {}'.format(pars["eMax"]))
        if pars["junctions"]:
            print('Minimum number of bins: {}'.format(pars["bins"]-5))
            print('Maximum number of bins: {}'.format(pars["bins"]+5))
            print('Junction points: {}'.format(pars["junctions"]))
        else:
            print('Number of bins: {}'.format(pars["bins"])) 
        print("****\n")

    
    binning = createBinning(opts, pars)
    '''
    if len(final_binning):
        with open(opts.output, "w") as out_binning:
            for elm in final_binning:
                out_binning.write(str(elm))
                out_binning.write("\n")
    else:
        print("No correct bin found. Nothing has been written to disk")
    '''

if __name__ == "__main__":
    main()
    print(datetime.now()-start)
