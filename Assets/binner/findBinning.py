#!/usr/bin/env python3
from datetime import datetime
import numpy as np

start = datetime.now()

def main(args=None):

    ## Binning parameters
    min_n_bins = 30
    max_n_bins = 50
    max_error = 2
    bins_vertex = [15, 100]

    for error in range(max_error,0,-1):
        best_error = [error, error]
        best_binning = []
        for bins in range(min_n_bins,max_n_bins):
            binning = np.logspace(0, 4, bins+1)
            found = [False, False]
            for bin_value in binning:
                if (np.absolute(bin_value-bins_vertex[0]))<=best_error[0]:
                    found[0] = True
                    if np.absolute(bin_value-bins_vertex[0]) < best_error[0]:
                        best_error[0] = np.absolute(bin_value-bins_vertex[0])
                if (np.absolute(bin_value-bins_vertex[1]))<=best_error[1]:
                    found[1] = True
                    if np.absolute(bin_value-bins_vertex[1]) < best_error[1]:
                        best_error[1] = np.absolute(bin_value-bins_vertex[1])
            if found[0]*found[1]==True:
                best_binning.append(binning)
        
        print('\nBinning found! Max error value: {} \tNumber of bins: {}'.format(error, len(best_binning[-1])-1 ))
        print("Best error: ", best_error)
        print(best_binning[-1])
    

if __name__ == "__main__":
    main()
    print(datetime.now()-start)
