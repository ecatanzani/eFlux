from argparse import ArgumentParser
import ROOT
from tqdm import tqdm

def GetROOTTreeName(filename: str) -> str:
    root_file = ROOT.TFile(filename,'READ')
    tree_name = None
    for key in root_file.GetListOfKeys():
        if root_file.Get(key.GetName()).ClassName() == 'TTree':
            tree_name = key.GetName()
            break
    return 

def Stats(input_list: str, verbose: bool) -> tuple:

    total_events = 0
    preselected_events = 0

    tree_name = None
    for file in tqdm(input_list, desc=f"Processing files in input directory: {input_list}"):
        if tree_name==None:
            tree_name = GetROOTTreeName(file)
        rdf = ROOT.RDataFrame(tree_name, file)
        total_events += rdf.Count().GetValue()
        preselected_events += rdf.Filter('evtfilter_all_cut==1').Count().GetValue()

    return (total_events, preselected_events)

def main(args=None):
    
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Flight DATA statistics")
    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='input directory')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)

    (total_events, preselected_events) = Stats(opts.input, opts.verbose)

    print(f"Total events: {total_events}")
    print(f"Preselected events: {preselected_events}")

if __name__ == '__main__':
    main()