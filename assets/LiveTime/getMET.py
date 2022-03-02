import ROOT
import datetime
from argparse import ArgumentParser

def get_date_from_data_file(_file):
    year = int(_file[_file.rfind('_OBS_')+5:_file.rfind('_OBS_')+9])
    month = int(_file[_file.rfind('_OBS_')+9:_file.rfind('_OBS_')+11])
    day = int(_file[_file.rfind('_OBS_')+11:_file.rfind('_OBS_')+13])

    return datetime.date(year, month, day)

def read_data_files(input_file_list, verbose):
    with open(input_file_list, "r") as _input_file:
        lines = _input_file.read().splitlines()
    if lines is not None:
        print('{} lines have been parsed from input file list [{}]'.format(len(lines), input_file_list))
        return lines, get_date_from_data_file(lines[0])

def extract_met(file_list, verbose):
    ROOT.gSystem.Load("libDmpService.so")
    evtch = ROOT.DmpChain("CollectionTree")
    for file in file_list:
        evtch.Add(file)
        if verbose:
            print('Adding {} to the chain...'.format(file))
    
    seconds = []
    for ev_idx in range(evtch.GetEntries()):
        pev = evtch.GetDmpEvent(ev_idx)
        evthdr = pev.pEvtHeader()
        seconds.append(int(evthdr.GetSecond()))
    
    return min(seconds), max(seconds)

def main(args=None):
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="Convert [Raw Data (2A)] from date to DAMPE MET")
    parser.add_argument("-i", "--input", type=str, dest='input', help='Input RAW data list')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False, action='store_true', help='run in high verbosity mode')
    
    opts = parser.parse_args(args)

    # Parse input data in data
    files, dampe_date = read_data_files(opts.input, opts.verbose)
    # Extract MET
    met_min, met_max = extract_met(files, opts.verbose)

    print('\n\n{} --> MET start: {} - MET end: {}\n'.format(dampe_date, met_min, met_max))
    


if __name__ == '__main__':
    main()