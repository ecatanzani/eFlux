import os
import subprocess
import ROOT
from argparse import ArgumentParser 

def findFiles(input_directory: str) -> list:
    files = []
    for folder in [os.path.join(input_directory, tmp_folder) for tmp_folder in os.listdir(input_directory) if tmp_folder.startswith('job_') and os.path.isdir(os.path.join(input_directory, tmp_folder))]:
        for file in [os.path.join(folder, 'outFiles', file) for file in os.listdir(os.path.join(folder, 'outFiles')) if file.endswith('.root')]:
            files.append(file)
    return files

def chainFiles(files: list, verbose: bool) -> ROOT.TChain:
    selection_chain = ROOT.TChain("electron_tree")
    for file in files:
        selection_chain.Add(file)
        if verbose:
            print(f"Adding {file} to the chain...")
    return selection_chain

def saveChain(selection_chain: ROOT.TChain, output_file: str, threads: int, verbose: bool):
    ROOT.ROOT.EnableImplicitMT(threads)
    selection_fr = ROOT.RDataFrame(selection_chain)
    selection_fr.Snapshot(selection_chain.GetName(), output_file)
    if verbose:
        print(f"Output file has been created: {output_file}")

def main(args=None):
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="Classified DAMPE DATA/MC selection extractor")
    parser.add_argument("-i", "--input_directory", type=str, 
                            dest='input_directory', help='input job directory', required=True)
    parser.add_argument("-o", "--output_file", type=str, 
                            dest='output_directory', help='output file', required=True)
    parser.add_argument("-v", "--verbose", dest='verbose', default=False, 
                            action='store_true', help='run in high verbosity mode')
    parser.add_argument("-t", "--threads", type=int, dest='threads',
                            const=1, nargs='?', help='Number of threads')
    opts = parser.parse_args(args)

    saveChain(chainFiles(findFiles(opts.input_directory), opts.verbose), opts.output_directory, opts.threads, opts.verbose)
    

if __name__ == '__main__':
    main()