#!/usr/bin/env python3
from datetime import datetime
from argparse import ArgumentParser
import numpy as np

start = datetime.now()


def createBinning(min_n_bins, max_n_bins, max_error, bins_vertex, opts):
    best_binning = []
    for error in range(max_error, 0, -1):
        best_error = []
        for idx in range(0, len(bins_vertex)):
            best_error.append(error)
        for bins in range(min_n_bins, max_n_bins):
            binning = np.logspace(0, 4, bins+1)
            found = [False, False]
            for bin_value in binning:
                if (np.absolute(bin_value-bins_vertex[0])) <= best_error[0]:
                    found[0] = True
                    if np.absolute(bin_value-bins_vertex[0]) < best_error[0]:
                        best_error[0] = np.absolute(bin_value-bins_vertex[0])
                if (np.absolute(bin_value-bins_vertex[1])) <= best_error[1]:
                    found[1] = True
                    if np.absolute(bin_value-bins_vertex[1]) < best_error[1]:
                        best_error[1] = np.absolute(bin_value-bins_vertex[1])
            if found[0]*found[1] == True:
                best_binning.append(binning)

    if opts.verbose:
        if best_binning:
            print("Binning found!")
            print("\n****")
            print("Max error value: ", error)
            print("Number of bins: ", len(best_binning[-1])-1)
            for idx, elm in enumerate(best_error):
                print('Best error bin vertex {}: {}'.format(idx, elm))
            print("****\n")
            print(best_binning[-1])
            print("\n")

    if best_binning:
        return best_binning[-1]
    else:
        return best_binning


def main(args=None):

    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Binning Finder")
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='Output binning text file')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)

    # Binning parameters
    min_n_bins = 30
    max_n_bins = 50
    max_error = 2
    bins_vertex = [15, 100]

    final_binning = createBinning(min_n_bins, max_n_bins, max_error, bins_vertex, opts)
    
    if len(final_binning):
        with open(opts.output, "w") as out_binning:
            for elm in final_binning:
                out_binning.write(str(elm))
                out_binning.write("\n")
    else:
        print("No correct bin found. Nothing has been written to disk")
    

if __name__ == "__main__":
    main()
    print(datetime.now()-start)
