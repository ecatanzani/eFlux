from argparse import ArgumentParser
import sys


def main(args=None):

    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Plot facility for Skimmed Flight Data")

    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output histo name')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')
    opts = parser.parse_args(args)
    
    # Load analysis functions
    sys.path.append("moduls")
    from configParser import parseConfigFile
    from downloadStatistics import getStatisticFiles
    from buildHisto import buildHisto
    from checkFiles import checkFiles

   
    # Get dictionary from config file parsing
    pars = parseConfigFile()
    # Download stats files
    files = getStatisticFiles(pars, opts)
    # Build final histo
    buildHisto(opts, files, pars)


if __name__ == "__main__":
    main()
