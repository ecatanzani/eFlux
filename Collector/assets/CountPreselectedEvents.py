from argparse import ArgumentParser
import ROOT
import os
from tqdm import tqdm

def GetROOTTreeName(filename: str) -> str:
    root_file = ROOT.TFile(filename,'READ')
    tree_name = None
    for key in root_file.GetListOfKeys():
        if root_file.Get(key.GetName()).ClassName() == 'TTree':
            tree_name = key.GetName()
            break
    return tree_name

def Stats(input_list: str, verbose: bool) -> tuple:

    with open(input_list) as file_list:
        input_files = file_list.read().splitlines()
    
    total_events = 0
    preselected_events = 0
    total_size = 0

    tree_name = None
    for file in tqdm(input_files, desc=f"Processing files in input list: {input_list}"):
        if tree_name==None:
            tree_name = GetROOTTreeName(file)
        rdf = ROOT.RDataFrame(tree_name, file)
        total_events += rdf.Count().GetValue()
        preselected_events += rdf.Filter('evtfilter_all_cut==1').Count().GetValue()
        total_size += os.path.getsize(file)

    return (total_events, preselected_events, total_size)

def main(args=None):
    
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Flight DATA statistics")
    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='input directory')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output file name')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)

    (total_events, preselected_events, total_size) = Stats(opts.input, opts.verbose)

    print(f"Total events: {total_events}")
    print(f"Preselected events: {preselected_events}")
    print(f"Total size (TB): {float(total_size)/1e-12}")

    with open(opts.output, 'w') as output:
        output.write(f"Total events: {total_events}\n")
        output.write(f"Preselected events: {preselected_events}\n")
        output.write(f"Preselection efficiency: {float(preselected_events)/total_events}\n")
        output.write(f"Total size (TB): {float(total_size)/1e-12}")
        
if __name__ == '__main__':
    main()