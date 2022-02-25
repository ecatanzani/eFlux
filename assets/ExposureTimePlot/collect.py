from argparse import ArgumentParser
from tqdm import tqdm
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os

def get_files(condor_input):
    csv_files = []
    for folder in os.listdir(condor_input):
        if folder.startswith('job_'):
            tmp_files = os.listdir(os.path.join(condor_input, folder, 'outFiles'))
            for file in tmp_files:
                if file.endswith('.csv'):
                    csv_files.append(os.path.join(condor_input, folder, 'outFiles', file))
    return csv_files

def buildPlot(infodb, output_file):
    
    # build the histo
    infodb["Date"]= pd.to_datetime(infodb["Date"])
    infodb_selected = infodb[(infodb['Date'] >= "2016-01-01") & (infodb['Date'] <= "2021-10-31")]
    fig, ax = plt.subplots()
    ax.plot(infodb_selected['Date'].values, infodb_selected['ExposureTime'].values, label="Exposure Time", color="cornflowerblue")
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%b"))
    
    _ = plt.xticks(rotation=30)

    fig.autofmt_xdate()
    
    plt.ylim(0, 8e+4)
    plt.xlabel('Date')
    plt.ylabel('Exposure Time (s)')

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
    complete_df = pd.concat([pd.read_csv(file) for file in files], ignore_index=True)
    complete_df = complete_df.sort_values(by="Date")
    complete_df.to_csv('{}.csv'.format(opts.output), index=False)
    # Build the plot
    if opts.output:
        buildPlot(complete_df, opts.output)

if __name__ == '__main__':
    main()