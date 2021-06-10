import os
from argparse import ArgumentParser

def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="DATA/MC list builder")
    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='input file directory')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output file list')
    
    opts = parser.parse_args(args)

    files = [f"{os.path.abspath(opts.input)}/{file}" for file in os.listdir(opts.input) if not file.startswith('.')]
    with open(opts.output, 'w') as outputlist:
        for file in files:
            outputlist.write(f"{file}\n")

if __name__ == '__main__':
    main()