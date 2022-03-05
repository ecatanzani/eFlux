#!/usr/bin/env python
import ROOT
from argparse import ArgumentParser

def compute_livetime(start_sec, end_sec):
    ROOT.gSystem.Load("libDmpService.so")
    return ROOT.DmpSvcLiveTime.GetInstance().GetLiveTime(start_sec, end_sec)

def main(args=None):
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="DAMPE live-time facility")
    parser.add_argument("-s", "--start", type=int, dest='start_met_time', help='start MET time')
    parser.add_argument("-e", "--end", type=int, dest='end_met_time', help='end MET time')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False, action='store_true', help='run in high verbosity mode')
    
    opts = parser.parse_args(args)

    if (opts.verbose):
        print('\nComputing LiveTime in the following interval: [{}, {}]\n'.format(opts.start_met_time, opts.end_met_time))
    
    print('LiveTime: {}'.format(compute_livetime(opts.start_met_time, opts.end_met_time)))


if __name__ == '__main__':
    main()