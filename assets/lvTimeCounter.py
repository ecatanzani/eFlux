#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser


def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="DAMPE live-time facility")
    parser.add_argument("-s", "--start", type=int, dest='start_time',
                        const=94608354, nargs='?', help='start time')
    parser.add_argument("-e", "--end", type=int, dest='end_time',
                        const=238981295, nargs='?', help='end time')
    opts = parser.parse_args(args)

    gSystem.Load("libDmpService.so")
    _time_instance = DmpSvcLiveTime.GetInstance()
    live_time = _time_instance.GetLiveTime(opts.start_time, opts.end_time)
    print('Live-Time: {}'.format(live_time))


if __name__ == '__main__':
    main()