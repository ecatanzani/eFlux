#!/usr/bin/env python
import time
import datetime
from ROOT import gSystem
from argparse import ArgumentParser

def data_to_timestamp(data):
    return int(time.mktime(datetime.datetime.strptime(data, "%d/%m/%Y").timetuple()))

def main(args=None):
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="DAMPE live-time facility")
    parser.add_argument("-s", "--start", type=str, dest='start_date', help='start date (dd/mm/yyyy)')
    parser.add_argument("-e", "--end", type=str, dest='end_date', help='start date (dd/mm/yyyy)')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False, action='store_true', help='run in high verbosity mode')
    
    opts = parser.parse_args(args)

    if opts.start_date is None:
        opts.start_date = "01/01/2016"
        if opts.verbose:
            print('The start date has not been provided... default one is used[{}]'.format(opts.start_date))
    if opts.end_date is None:
        opts.end_date = "03/11/2021"
        if opts.verbose:
            print('The end date has not been provided... default one is used[{}]'.format(opts.end_date))

    if (opts.verbose):
        print('Computing LiveTime in the following interval: [{}, {}]'.format(opts.start_date, opts.end_date))
    gSystem.Load("libDmpService.so")
    _time_instance = DmpSvcLiveTime.GetInstance(data_to_timestamp(opts.start_date), data_to_timestamp(opts.end_date))
    live_time = _time_instance.GetLiveTime()
    
    print('LiveTime: {}'.format(live_time))


if __name__ == '__main__':
    main()