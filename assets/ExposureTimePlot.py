import os
import ROOT
from tqdm import tqdm
import pandas as pd
import numpy as np
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib import rcParams
from datetime import date

def GetROOTTreeName(filename):
    root_file = ROOT.TFile(filename,'READ')
    tree_name = None
    for key in root_file.GetListOfKeys():
        if root_file.Get(key.GetName()).ClassName() == 'TTree':
            tree_name = key.GetName()
            break
    return tree_name

def extract_seconds_from_file(file, tree_name, verbose):
    evtch = ROOT.DmpChain(tree_name)
    evtch.Add(file)
    if verbose:
        print('Adding {} to the chain...'.format(file))
    
    seconds = np.empty(evtch.GetEntries(), dtype=int)
    for ev_idx in range(evtch.GetEntries()):
        pev = evtch.GetDmpEvent(ev_idx)
        evthdr = pev.pEvtHeader()
        seconds[ev_idx] = int(evthdr.GetSecond())
    
    return np.min(seconds), np.max(seconds)

def compute_exposure(start_sec, end_sec):
    return ROOT.DmpSvcLiveTime.GetInstance().GetLiveTime(start_sec, end_sec)

def ParseInputFiles(input_dir, verbose):
    
    df = pd.DataFrame(columns=['Date', 'ExposureTime'])

    # Sort folders in order to have an ordinate dictionary of values
    folders = os.listdir(input_dir)
    folders.sort()

    for folder in tqdm(folders, desc='Processing folders in {}'.format(input_dir)):
        
        year = int(folder[:4])
        month = int(folder[4:6])
        day = int(folder[6:8])
        
        root_tree_name = None

        if folder.startswith('20') and os.path.isfile(os.path.join(input_dir, folder, 'output.log')):
            
            # Load the RDF
            root_filename = [os.path.join(input_dir, folder, "outFiles", file) for file in os.listdir(os.path.join(input_dir, folder, "outFiles")) if file.endswith('.root')][0]
            if root_tree_name == None:
                root_tree_name = GetROOTTreeName(root_filename)

            # Create new row in DataFrame
            df = df.append(pd.Series({'Date': date(year, month, day), 'ExposureTime': compute_exposure(extract_seconds_from_file(root_filename, root_tree_name, verbose))}))

    return df
    
def buildPlot(infodb, output_file):
    # build the histo
    rcParams.update({'figure.autolayout': True})
    
    plt.plot(infodb['Date'].values, infodb['ExposureTime'].values, label="Exposure Time", color="cornflowerblue")

    plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
    plt.yscale('log')
    plt.ylim(10, 1e+8)
    plt.ylabel("counts/day")
    plt.savefig('{}.pdf'.format(output_file))

def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Flight DATA statistics")
    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='input directory')
    parser.add_argument("-c", "--input-csv", type=str,
                        dest='input_csv', help='input csv file')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output file name')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)
    
    info = pd.DataFrame()

    if opts.input:
        # Load DAMPE Service Library
        ROOT.gSystem.Load("libDmpService.so")
        
        # Parse input files per day and compute exposure time
        info = ParseInputFiles(opts.input, opts.verbose)
        info.to_csv('{}.csv'.format(opts.output))
    if opts.input_csv:
        info = pd.read_csv(opts.input_csv)
    
    # Build the plot
    buildPlot(info, opts.output)

if __name__ == '__main__':
    main()