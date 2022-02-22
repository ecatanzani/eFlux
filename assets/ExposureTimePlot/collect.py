from argparse import ArgumentParser
from tqdm import tqdm
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
import os

def get_files(condor_input):
    csv_files = []
    for folder in os.listdir(condor_input) if folder.startswith('job_'):
        tmp_files = os.listdir(os.path.join(condor_input, folder, 'outFiles'))
        for file in tmp_files if file.endswith('.csv'):
            csv_files.append(os.path.join(condor_input, folder, 'outFiles', file))
    return csv_files

def buildPlot(infodb, output_file):
    register_matplotlib_converters()
    # build the histo
    matplotlib.rcParams.update({'figure.autolayout': True})
    plt.plot(infodb['Date'].values, infodb['ExposureTime'].values, label="Exposure Time", color="cornflowerblue")
    plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
    plt.yscale('log')
    plt.ylim(1e+3, 1e+6)
    plt.ylabel("counts/day")
    plt.savefig('{}.pdf'.format(output_file))

def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Flight DATA statistics - Condor Engine")
    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='input HTCondor directory')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output file')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)
    
    # Read output csv files
    files = get_files(opts.input)
    # Build final dataframe and sort by date
    complete_df = pd.concat(files, ignore_index=True)
    complete_df = complete_df.sort_values(by="Date")
    # Build the plot
    if opts.output:
        buildPlot(complete_df, opts.output)

if __name__ == '__main__':
    main()