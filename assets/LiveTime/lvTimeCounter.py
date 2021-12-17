#!/usr/bin/env python
import time
import datetime
import ROOT
from argparse import ArgumentParser

def data_to_timestamp(data):
    return int(time.mktime(datetime.datetime.strptime(data, "%d/%m/%Y").timetuple()))

def compute_livetime(start_sec, end_sec):
    return ROOT.DmpSvcLiveTime.GetInstance().GetLiveTime(start_sec, end_sec)

def extract_seconds_from_files(files, verbose):
    evtch = ROOT.DmpChain("CollectionTree")
    for single_file in files:
        evtch.Add(single_file)
        if verbose:
            print('Adding {} to the chain...'.format(single_file))
    
    seconds = []
    for ev_idx in range(evtch.GetEntries()):
        pev = evtch.GetDmpEvent(ev_idx)
        evthdr = pev.pEvtHeader()
        seconds.append(int(evthdr.GetSecond()))
    
    return min(seconds), max(seconds)

def main(args=None):
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="DAMPE live-time facility")
    parser.add_argument("-s", "--start", type=str, dest='start_file', help='start 2A DAMPE data file')
    parser.add_argument("-e", "--end", type=str, dest='end_file', help='end 2A DAMPE data file')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False, action='store_true', help='run in high verbosity mode')
    
    opts = parser.parse_args(args)
    ROOT.gSystem.Load("libDmpService.so")
    start_second, end_second = extract_seconds_from_files((opts.start_file, opts.end_file), opts.verbose)

    if (opts.verbose):
        print('\nComputing LiveTime in the following interval: [{}, {}]\n'.format(start_second, end_second))
    
    print('LiveTime: {}'.format(compute_livetime(start_second, end_second)))


if __name__ == '__main__':
    main()