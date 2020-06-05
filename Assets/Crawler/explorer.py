#!/usr/bin/env python3
from datetime import datetime
from argparse import ArgumentParser
import sys

start = datetime.now()


def main(args=None):

    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="reco MC files explorer")

    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='json data storage')
    parser.add_argument("-f", "--farm", type=str,
                        dest='farm', help='MC ROOT input farm')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output file list')
    parser.add_argument("-d", "--data", dest='data', default=False,
                        action='store_true', help='get DATA file list')
    parser.add_argument("-m", "--mc", dest='mc', default=False,
                        action='store_true', help='get MC file list')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')
    parser.add_argument("-s", "--spectrum", type=str,
                        dest='spectrum', help='Study of the energy spectrum of a given dataSet')

    opts = parser.parse_args(args)

    # Load analysis functions
    sys.path.append("moduls")
    from configParser import parseConfigFile
    from DATA_list import createDATAlist
    from MC_list import createMClist
    from singleSample import listMCsample

    # Get dictionary from config file parsing
    pars = parseConfigFile(opts)

    if opts.data:
        createDATAlist(pars, opts)
    if opts.mc:
        createMClist(pars, opts)
    if opts.spectrum:
        listMCsample(pars, opts)


if __name__ == "__main__":
    main()
    print(datetime.now()-start)
