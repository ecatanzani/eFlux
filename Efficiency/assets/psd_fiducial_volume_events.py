import os
import argparse
from tqdm import tqdm
from ROOT import TFile

def TestROOTFile(path: str) -> bool:
    _tmp_file = TFile(path)
    if _tmp_file and not _tmp_file.IsOpen():
        return False
    elif _tmp_file and _tmp_file.IsOpen() and _tmp_file.IsZombie():
        _tmp_file.Close()
        return False
    elif _tmp_file and _tmp_file.IsOpen() and _tmp_file.TestBit(TFile.kRecovered):
        _tmp_file.Close()
        return False
    else:
        _tmp_file.Close()
        return True

def getListOfFiles(condor_wd: str) -> list:
    data_dirs = []
    for tmp_dir in tqdm(os.listdir(condor_wd), desc='Scanning local HTCondor dir'):
        if tmp_dir.startswith('job_'):
            full_dir_path = os.path.join(condor_wd, tmp_dir)
            expected_condor_outDir = os.path.join(full_dir_path, "outFiles")
            # Check if 'outFiles' dir exists
            if os.path.isdir(expected_condor_outDir):
                _list_dir = [f"{expected_condor_outDir}/{file}" for file in os.listdir(expected_condor_outDir) if file.endswith('.root')]
                skipped_dir = False
                for file_idx, tmp_acc_full_path in enumerate(_list_dir):
                    # Check if output ROOT file exists
                    if os.path.isfile(tmp_acc_full_path):
                        # Check if output ROOT file is redable
                        if TestROOTFile(tmp_acc_full_path):
                            tmp_acc_file = TFile.Open(tmp_acc_full_path, "READ")
                            # Check if output ROOT file is redable
                            if tmp_acc_file.IsOpen():
                                # Check if output ROOT file has keys
                                outKeys = tmp_acc_file.GetNkeys()
                                if outKeys:
                                    if file_idx == len(_list_dir)-1 and not skipped_dir:
                                        data_dirs.append(full_dir_path)
                                else:
                                    # output ROOT file has been open but has not keys
                                    if not skipped_dir:
                                        skipped_dir = True
                        else:
                            # output ROOT file has not been opened correctly
                            if not skipped_dir:
                                skipped_dir = True
                    else:
                        # output ROOT file does not exist
                        if not skipped_dir:
                            skipped_dir = True
    
    return [os.path.join(file, "output.log") for file in data_dirs]

def get_log_stats(log_files: list) -> tuple:
    bgo_track_preselection = 0
    bgo_track_preselection_psd_fvolume = 0
    bgo_track_preselection_no_psd_fvolume = 0

    for log_file in log_files:
        with open(log_file, 'r') as f:
            for line in f:
                if 'Number of events after BGO and STK selection:' in line:
                    bgo_track_preselection += int(line[line.rfind(': ')+2:])
                if 'Number of events within the PSD fiducial volume:' in line:
                    bgo_track_preselection_psd_fvolume += int(line[line.rfind(': ')+2:])
                if 'Number of events outside the PSD fiducial volume:' in line:
                    bgo_track_preselection_no_psd_fvolume += int(line[line.rfind(': ')+2:])

    return (bgo_track_preselection, bgo_track_preselection_psd_fvolume, bgo_track_preselection_no_psd_fvolume)

def main(args=None):
    parser = argparse.ArgumentParser(description='PSD Study Facility')

    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='Input DATA/MC folder')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output text file')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)

    if opts.input:
        log_files = getListOfFiles(opts.input)
        stats = get_log_stats(log_files)

        print(f"Number of events after BGO and STK selection: {stats[0]}")
        print(f"Number of events within the PSD fiducial volume: {stats[1]}")
        print(f"Number of events outside the PSD fiducial volume: {stats[2]}")

        if opts.output:
            with open(opts.output, 'w') as f:
                f.write(f"Number of events after BGO and STK selection: {stats[0]}\n")
                f.write(f"Number of events within the PSD fiducial volume: {stats[1]}\n")
                f.write(f"Number of events outside the PSD fiducial volume: {stats[2]}\n")

if __name__ == '__main__':
    main()