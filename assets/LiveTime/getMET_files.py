import ROOT
import datetime
import subprocess
from argparse import ArgumentParser

def get_date(input_date: str, verbose: bool) -> datetime.date:
    year = int(input_date[0:4])
    month = int(input_date[4:6])
    day = int(input_date[6:8])
    
    if verbose:
        print(f"Searching RAW data files for {year}-{month}-{day}")
    
    return datetime.date(year, month, day)

def get_date_from_dir(_dir: str) -> datetime.date:
    year = int(_dir[_dir.rfind('/')+1:_dir.rfind('/')+5])
    month = int(_dir[_dir.rfind('/')+5:_dir.rfind('/')+7])
    day = int(_dir[_dir.rfind('/')+7:_dir.rfind('/')+9])
    
    return datetime.date(year, month, day)

def get_date_from_data_file(_file: str) -> datetime.date:
    year = int(_file[_file.rfind('_OBS_')+5:_file.rfind('_OBS_')+9])
    month = int(_file[_file.rfind('_OBS_')+9:_file.rfind('_OBS_')+11])
    day = int(_file[_file.rfind('_OBS_')+11:_file.rfind('_OBS_')+13])

    return datetime.date(year, month, day)

def check_date(_dir: str, input_date: datetime.date) -> bool:
    start_date = input_date - datetime.timedelta(days=1)
    end_date = input_date + datetime.timedelta(days=1)
    tmp_date = get_date_from_dir(_dir)

    return True if start_date <= tmp_date <= end_date else False

def get_data_files(input_date: datetime.date, verbose: bool) -> list:
    files = []
    pars = {
        'farmAddress': 'root://xrootd-dampe.cloud.ba.infn.it//',
        'data_XRDFS_path': '/FM/FlightData/2A/'
    }

    # Get stage 0 dirs --> /FM/FlightData/2A/
    _cmd = 'xrdfs {} ls {}'.format(pars['farmAddress'], pars['data_XRDFS_path'])
    _out = subprocess.run(_cmd, shell=True, check=True, stdout=subprocess.PIPE)
    _dirs = str.split(_out.stdout.decode('utf-8').rstrip(), '\n')

    # Get stage 1 dirs --> /FM/FlightData/2A/DayOfAcquisition/
    for _dir in _dirs:
        if '/2A/20' in _dir:
            if check_date(_dir, input_date):
                
                if verbose:
                    print(f"Scanning [{_dir}] ...")

                _cmd = 'xrdfs {} ls {}'.format(pars['farmAddress'], _dir)
                _out = subprocess.run(_cmd, shell=True, check=True, stdout=subprocess.PIPE)
                _acquisition_dirs = str.split(_out.stdout.decode('utf-8').rstrip(), '\n')

                # Get stage 2 dirs --> /FM/FlightData/2A/DayOfAcquisition/DataDir
                for _acquisition_dir in _acquisition_dirs:
                    _cmd = 'xrdfs {} ls {}'.format(pars['farmAddress'], _acquisition_dir)
                    _out = subprocess.run(_cmd, shell=True, check=True, stdout=subprocess.PIPE)
                    _acquisition_files = str.split(_out.stdout.decode('utf-8').rstrip(), '\n')

                    # Get ROOT data files
                    for _acquisition_file in _acquisition_files:
                        if _acquisition_file.endswith('.root'):
                            files.append(f"{pars['farmAddress']}{_acquisition_file}")

    return files

def filter_data_files(file_list: list, input_date: datetime.date, verbose: bool) -> list:
    filtered_files = []
    for file in file_list:
        if get_date_from_data_file(file) == input_date:
            filtered_files.append(file)
    
    if verbose:
        print(f"\nFound {len(filtered_files)} files over {len(file_list)} for {input_date}\n")
    return filtered_files

def main(args=None):
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="Convert [Raw Data (2A)] from date to DAMPE MET")
    parser.add_argument("-d", "--date", type=str, dest='date', help='date: yyyy-mm-dd')
    parser.add_argument("-o", "--output", type=str, dest='output', help='output_file')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False, action='store_true', help='run in high verbosity mode')
    
    opts = parser.parse_args(args)

    # Parse input data in datetime format
    dampe_date = get_date(opts.date, opts.verbose)
    # Extract RAW data files considering the input date +- 1 day 
    file_list = get_data_files(dampe_date, opts.verbose)
    # Parse only the files belonging to the input date
    filtered_files = filter_data_files(file_list, dampe_date, opts.verbose)
    if (opts.verbose):
        print("\n*** List of filtered RAW data files ***\n")
        for file in filtered_files:
            print(file)
    
    with open(opts.output, 'w') as _output_file:
        for file in filtered_files:
            _output_file.write(f"{file}\n")

if __name__ == '__main__':
    main()