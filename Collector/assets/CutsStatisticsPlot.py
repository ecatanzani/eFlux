import os
import ROOT
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib import rcParams
from datetime import date, timedelta

def GetROOTTreeName(filename: str) -> str:
    root_file = ROOT.TFile(filename,'READ')
    tree_name = None
    for key in root_file.GetListOfKeys():
        if root_file.Get(key.GetName()).ClassName() == 'TTree':
            tree_name = key.GetName()
            break
    return tree_name


def BuildPlot(input_dir: str, output_file: str, verbose: bool):
    stats = {
        'dates': [],
        'raw_events': [],
        'bgo_fiducial': [],
        'nbarlayer13': [],
        'maxrms': [],
        'track_selection': [],
        'psd_stk_match': [],
        'psd_charge': []
    }

    for folder in os.listdir(input_dir):
        
        year = int(folder[:4])
        month = int(folder[4:6])
        day = int(folder[6:8])
        
        if folder.startswith('20'):

            # Extract the date information from the folder name
            stats['dates'].append(date(year, month, day))

            # Load the RDF
            root_filename = [file for file in os.listdir(f"{input_dir}/{folder}/outFiles/") if file.endswith('.root')][0]
            root_tree_name = GetROOTTreeName(root_filename)

            # Extract the counts info
            if root_tree_name is not None:
                rdf = ROOT.RDataFrame(root_tree_name, root_filename)
                stats['raw_events'].append(rdf.Count().GetValue())
                stats['bgo_fiducial'].append(rdf.Filter('evtfilter_BGO_fiducial==1').Count().GetValue())
                stats['nbarlayer13'].append(rdf.Filter('evtfilter_nBarLayer13_cut==1').Count().GetValue())
                stats['maxrms'].append(rdf.Filter('evtfilter_maxRms_cut==1').Count().GetValue())
                stats['track_selection'].append(rdf.Filter('evtfilter_track_selection_cut==1').Count().GetValue())
                stats['psd_stk_match'].append(rdf.Filter('evtfilter_psd_stk_match_cut==1').Count().GetValue())
                stats['psd_charge'].append(rdf.Filter('evtfilter_psd_charge_cut==1').Count().GetValue())

    # build the histo
    rcParams.update({'figure.autolayout': True})
    plt.plot(stats['dates'], stats['raw_events'], label="total counts", color="dimgray")
    plt.plot(stats['dates'], stats['bgo_fiducial'], label="BGO fiducial cut", color="cornflowerblue")
    plt.plot(stats['dates'], stats['nbarlayer13'], label="nBarLayer13 cut", color="darkorange")
    plt.plot(stats['dates'], stats['maxrms'], label="max RMS cut", color="forestgreen")
    plt.plot(stats['dates'], stats['track_selection'], label="Track Selection cut", color="crimson")
    plt.plot(stats['dates'], stats['psd_stk_match'], label="PSD/STK match cut", color="blueviolet")
    plt.plot(stats['dates'], stats['psd_charge'], label="PSD charge cut", color="saddlebrown")

    plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
    plt.yscale('log')
    plt.ylim(10, 1e+8)
    plt.ylabel("counts/day")
    plt.savefig(output_file)
    

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
    
    BuildPlot(opts.input, opts.output, opts.verbose)

if __name__ == '__main__':
    main()