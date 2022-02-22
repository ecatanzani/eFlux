import ROOT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
import numpy as np
import pandas as pd
from tqdm import tqdm
from datetime import date
from argparse import ArgumentParser

def ExtractDateFromFile(file):
    year = int(file[file.rfind('_OBS_')+5:file.rfind('_OBS_')+9])
    month = int(file[file.rfind('_OBS_')+9:file.rfind('_OBS_')+11])
    day = int(file[file.rfind('_OBS_')+11:file.rfind('_OBS_')+13])
    return date(year, month, day)

def ParseInputFiles(input_files):
    with open(input_files, 'r') as inputlist:
        filelist = inputlist.read().splitlines()
    date = None
    day_file_list = []
    day_info = {}
    for file in tqdm(filelist, desc='Parsing input list...'):
        tmp_date =  ExtractDateFromFile(file)
        if date is None:
            date = tmp_date
        else:
            if date != tmp_date:
                day_info[date] = day_file_list
                date = tmp_date
                day_file_list = []
        day_file_list.append(file)
    if day_file_list:
        day_info[date] = day_file_list
    return day_info       

def extract_seconds_from_files(file_list, verbose, file_tree_name='CollectionTree'):
    if verbose:
        print('Extracting seconds from files...')
    evtch = ROOT.DmpChain(file_tree_name)
    for file in file_list:
        evtch.Add(file)
        if verbose:
            print('Adding {} to the chain...'.format(file))
    seconds = np.empty(evtch.GetEntries(), dtype=int)
    for ev_idx in range(evtch.GetEntries()):
        pev = evtch.GetDmpEvent(ev_idx)
        evthdr = pev.pEvtHeader()
        seconds[ev_idx] = int(evthdr.GetSecond())
    min_sec = np.min(seconds)
    max_sec = np.max(seconds)
    if verbose:
        print('Extracted seconds from files: (min) {} - (max) {}'.format(min_sec, max_sec))
    return min_sec, max_sec

def GetTimeEdges(day_info, verbose):
    max_min_sec = []
    for tmp_date in day_info:
        min_sec, max_sec = extract_seconds_from_files(day_info[tmp_date], verbose)
        max_min_sec.append((min_sec, max_sec))
    return max_min_sec

def BuildDataFrame(day_info, verbose):
    max_min_sec = GetTimeEdges(day_info, verbose)
    df = pd.DataFrame(columns=['Date', 'ExposureTime'])
    for idx, tmp_date in tqdm(enumerate(day_info), desc='Building dataframe...'):
        exposure = ROOT.DmpSvcLiveTime.GetInstance().GetLiveTime(max_min_sec[idx][0], max_min_sec[idx][1])
        df = df.append(pd.Series({'Date': tmp_date, 'ExposureTime': exposure}), ignore_index=True)
    # Sort DataFrame by date
    df = df.sort_values(by="Date")
    return df
    
def buildPlot(infodb, output_file):
    register_matplotlib_converters()
    # build the histo
    matplotlib.rcParams.update({'figure.autolayout': True})
    plt.plot(infodb['Date'].values, infodb['ExposureTime'].values, label="Exposure Time", color="cornflowerblue")
    plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
    plt.yscale('log')
    plt.ylim(1e+3, 1e+6)
    plt.ylabel("counts/day")
    plt.savefig('{}.pdf'.format(output_file))

def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Flight DATA statistics")
    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='input DATA list')
    parser.add_argument("-c", "--input-csv", type=str,
                        dest='input_csv', help='input csv file')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output file name')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)
    
    info = pd.DataFrame()

    # Load DAMPE Service Library
    ROOT.gSystem.Load("libDmpService.so")

    if opts.input:
        # Parse input files per day and compute exposure time
        info = BuildDataFrame(ParseInputFiles(opts.input), opts.verbose)
        info.to_csv('{}.csv'.format(opts.output))
    if opts.input_csv:
        info = pd.read_csv(opts.input_csv)
    
    # Build the plot
    buildPlot(info, opts.output)

if __name__ == '__main__':
    main()