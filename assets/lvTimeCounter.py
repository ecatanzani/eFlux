#!/usr/bin/env python
import time
import datetime
from ROOT import *
from argparse import ArgumentParser

def data_to_timestamp(data):
    return int(time.mktime(datetime.datetime.strptime(data, "%d/%m/%Y").timetuple()))

def compute_livetime(start_date, end_date):
    gSystem.Load("libDmpService.so")
    return DmpSvcLiveTime.GetInstance().GetLiveTime(data_to_timestamp(start_date), data_to_timestamp(end_date))

def main(args=None):
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="DAMPE live-time facility")
    parser.add_argument("-s", "--start", type=str, dest='start_date', help='start date (dd/mm/yyyy)')
    parser.add_argument("-e", "--end", type=str, dest='end_date', help='start date (dd/mm/yyyy)')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False, action='store_true', help='run in high verbosity mode')
    
    opts = parser.parse_args(args)

    if opts.start_date is None:
        opts.start_date = "01/01/2016"
        if opts.verbose:
            print('\nThe start date has not been provided... default one is used[{}]'.format(opts.start_date))
    if opts.end_date is None:
        opts.end_date = "03/11/2021"
        if opts.verbose:
            print('\nThe end date has not been provided... default one is used[{}]'.format(opts.end_date))

    if (opts.verbose):
        print('\n\nComputing LiveTime in the following interval: [{}, {}]\n\n'.format(opts.start_date, opts.end_date))
    
    print('LiveTime: {}'.format(compute_livetime(opts.start_date, opts.end_date)))


if __name__ == '__main__':
    main()