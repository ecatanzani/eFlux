import os
import subprocess
from argparse import ArgumentParser

def findBinFiles(input_directory: str, n_bins: int) -> list:
    files = []
    for folder in [os.path.join(input_directory, tmp_folder) for tmp_folder in os.listdir(input_directory) if tmp_folder.startswith('job_') and os.path.isdir(os.path.join(input_directory, tmp_folder))]:
        for file in [file for file in os.listdir(os.path.join(folder, 'outFiles')) if '_energybin_' in file and file.endswith('.root')]:
            file.append(file)
    
    split_bins = []
    for bin in range(1, n_bins+1):
        bin_files = [file for file in files if f"_energybin_{bin}" in file]
        split_bins.append(bin_files)

    return split_bins

def aggregate(data_files: list, input_dir: str, output_dir: str, add_file_number: int, verbose: bool):
    if verbose:
        print(f"Going to add {len(data_files)} ROOT files...")
    
    _file_list = str()
    _list_idx = 0
    _tmp_out_name = os.path.join(output_dir, input_dir[input_dir.rfind('/')+1:])

    for fidx, file in enumerate(data_files):
        _file_list += f" {file}"
        if (fidx+1) % add_file_number == 0:
            _out_full_name = f"{_tmp_out_name}_{_list_idx}.root"
            _cmd = f"hadd {_out_full_name}{_file_list}"
            if verbose:
                print(_cmd)
            subprocess.run(_cmd, shell=True, check=True)
            _file_list = f" {_out_full_name}"
            _list_idx += 1
    if _file_list:
        _out_full_name = f"{_tmp_out_name}_{_list_idx}.root"
        _cmd = f"hadd {_out_full_name}{_file_list}"
        if verbose:
            print(_cmd)
        subprocess.run(_cmd, shell=True, check=True)

def addFiles(file_list: list, input_dir: str, output_dir: str, add_number: int, verbose: bool):
    for files in file_list:
        aggregate(files, input_dir, output_dir, add_number, verbose)


def main(args=None):
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="Classified DAMPE DATA integrator")
    parser.add_argument("-i", "--input_directory", type=str, 
                            dest='input_directory', help='input job directory', required=True)
    parser.add_argument("-o", "--output_directory", type=str, 
                            dest='output_directory', help='output job directory', required=True)
    parser.add_argument("-v", "--verbose", dest='verbose', default=False, 
                            action='store_true', help='run in high verbosity mode')
    parser.add_argument("-b", "--bins", type=int, dest='bins',
                            const=50, nargs='?', help='Number of bins')
    parser.add_argument("-n", "--number", type=int, dest='number',
                            const=10, nargs='?', help='Number of files to add')
    opts = parser.parse_args(args)

    addFiles(findBinFiles(opts.input_directory, opts.bins), opts.input_directory, opts.output_directory, opts.number, opts.verbose)
    

if __name__ == '__main__':
    main()