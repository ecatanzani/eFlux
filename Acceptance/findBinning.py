#!/usr/bin/env python3

import numpy as np

input = {
    'min_energy' : 1,
    'max_energy' : 20000,
    'num' : 50,
    'steps' : [10, 100, 1000, 10000],
    'tolerance' : 30
}

def isBinningAcceptable(bins: np.ndarray) -> bool:
    distances = []
    good_binning = True
    for elm in input['steps']:
        min_distance = min([abs(elm-tmp_elm) for tmp_elm in bins])
        if min_distance > input['tolerance']:
            good_binning = False
            break
        else:
            distances.append(min_distance)
    return good_binning

def inputBinning() -> bool:
    return isBinningAcceptable(np.logspace(np.log10(input['min_energy']), np.log10(input['max_energy']), input['num']))

def main(args=None):

    if inputBinning():
        print(f"\nGood binning has been found with parameters: {input}")
        print(f"\n{np.logspace(np.log10(input['min_energy']), np.log10(input['max_energy']), input['num'])}")
    else:
        good_binning = False
        num = input['num']-10
        while good_binning==False:
            if num>input['num']+10:
                break
            bins = np.logspace(np.log10(input['min_energy']), np.log10(input['max_energy']), num)
            good_binning = isBinningAcceptable(bins)
            if not good_binning:
                num += 1
        if good_binning:
            input['num'] = num
            print(f"\nGood binning has been found with parameters: {input}")
            print(f"\n{np.logspace(np.log10(input['min_energy']), np.log10(input['max_energy']), input['num'])}")
        else:
            print('No succesful match has been found... try with different parameters')
        
if __name__ == '__main__':
    main()