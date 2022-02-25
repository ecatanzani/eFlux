from argparse import ArgumentParser
import subprocess
import shutil
from tqdm import tqdm
import pandas as pd
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

def check_db(files):

    for file in files:
        df = pd.read_csv(file)
        if 0 in df['ExposureTime'].values:
            print("{} has no exposure time ... resubmitting ...".format(file))

            output_files = [os.path.join(file[:file.rfind('/outFiles/')], rem_file) for rem_file in os.listdir(file[:file.rfind('/outFiles/')]) if rem_file.startswith('output')]
            outfiles_dir = os.path.join(file[:file.rfind('/outFiles/')], 'outFiles')

            for output_file in output_files:
                try:
                    os.remove(output_file)
                except OSError as e:
                    print("Error: %s : %s" % (output_file, e.strerror))
            
            if os.path.exists(outfiles_dir):
                try:
                    shutil.rmtree(outfiles_dir)
                except OSError as e:
                    print("Error: %s : %s" % (outfiles_dir, e.strerror))

            condor_subfile = os.path.join(file[:file.rfind('/outFiles/')], 'cndr.sub')
            if os.path.exists(condor_subfile):
                print(condor_subfile)
                subprocess.run(['condor_submit -name sn-02.cr.cnaf.infn.it -spool {}'.format(condor_subfile)], shell=True, check=True)

def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Flight DATA statistics - Condor Engine")
    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='input HTCondor directory')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)
    
    # Read output csv files
    files = get_files(opts.input)
    # Build final dataframe and sort by date
    check_db(files)

if __name__ == '__main__':
    main()