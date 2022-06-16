import os
import datetime
import subprocess
from tqdm import tqdm
from argparse import ArgumentParser

def ExtractDateFromFile(file: str) -> datetime.date:
    '''
    year = int(file[file.rfind('/2A/')+4:file.rfind('/2A/')+8])
    month = int(file[file.rfind('/2A/')+8:file.rfind('/2A/')+10])
    day = int(file[file.rfind('/2A/')+10:file.rfind('/2A/')+12])
    '''
    year = int(file[file.rfind('_OBS_')+5:file.rfind('_OBS_')+9])
    month = int(file[file.rfind('_OBS_')+9:file.rfind('_OBS_')+11])
    day = int(file[file.rfind('_OBS_')+11:file.rfind('_OBS_')+13])
    return datetime.date(year, month, day)

def WiteListOnFile(list: list, list_date: datetime.date, output: str) -> str:

    month = str(list_date.month) if list_date.month > 9 else f"0{list_date.month}"
    day = str(list_date.day) if list_date.day > 9 else f"0{list_date.day}"
    newdir = f"{output}/{list_date.year}{month}{day}"
    os.mkdir(newdir)
    with open(f"{newdir}/dataDayList.txt", 'w') as outputlist:
        for file in list:
            outputlist.write(f"{file}\n")
    return newdir

def ParseList(pars: dict) -> list:

    with open(pars['input_data_list'], 'r') as inputlist:
        filelist = inputlist.read().splitlines()
    
    date = None
    day_file_list = []
    day_folders = []
    for file in tqdm(filelist, desc='Parsing input list...'):
        tmp_date =  ExtractDateFromFile(file)
        if date is None:
            date = tmp_date
        else:
            if date != tmp_date:
                day_folders.append(WiteListOnFile(day_file_list, date, pars['output']))
                date = tmp_date
                day_file_list = []
        day_file_list.append(file)
    
    if day_file_list:
        day_folders.append(WiteListOnFile(day_file_list, date, pars['output']))
    
    return day_folders

def WriteCondorFiles(day_folders: list):
    for folder in day_folders:

        outputPath = f"{folder}/output.log"
        logPath = f"{folder}/output.clog"
        errPath = f"{folder}/output.err"
        bashScriptPath = f"{folder}/script.sh"
        subFilePath = f"{folder}/cndr.sub"

        with open(subFilePath, 'w') as condor_submit_file:
            condor_submit_file.write("universe = vanilla\n")
            condor_submit_file.write(f"executable = {bashScriptPath}\n")
            condor_submit_file.write(f"output = {outputPath}\n")
            condor_submit_file.write(f"error = {errPath}\n")
            condor_submit_file.write(f"log = {logPath}\n")
            condor_submit_file.write("ShouldTransferFiles = YES\n")
            condor_submit_file.write("WhenToTransferOutput = ON_EXIT\n")
            condor_submit_file.write("queue 1")

def WriteCondorScripts(day_folders: list, pars: dict):
    for folder in day_folders:
        
        dataListPath = f"{folder}/dataDayList.txt"
        bashScriptPath = f"{folder}/script.sh"
        tmpOutDir = f"{folder}/outFiles"
        
        with open(bashScriptPath, 'w') as bash_script_file:
            bash_script_file.write("#!/usr/bin/env bash\n")
            bash_script_file.write("source /cvmfs/dampe.cern.ch/centos7/etc/setup.sh\n")
            bash_script_file.write("dampe_init trunk\n")
            bash_script_file.write(f"mkdir {tmpOutDir}\n")
            _cmd = f"{pars['executable']} -w {pars['config']} -i {dataListPath} -d {tmpOutDir} -v -r"
            bash_script_file.write(_cmd)

def SubmitJobs(day_folders: list):
    for folder in day_folders:
        subFilePath = f"{folder}/cndr.sub"
        subprocess.run([f"condor_submit -name sn-02.cr.cnaf.infn.it -spool {subFilePath}"], shell=True, check=True)

def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Flight DATA statistics")
    parser.add_argument("-l", "--list", type=str,
                        dest='list', help='input file list')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output directory')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    parser.add_argument("-c", "--config", type=str,
                        dest='config', help='Software Config Directory')
    parser.add_argument("-x", "--executable", type=str,
                        dest='executable', help='Analysis script')

    opts = parser.parse_args(args)

    pars = {
        "input_data_list": opts.list,
        "config": opts.config,
        "output": opts.output,
        "executable": opts.executable,
        "verbose": opts.verbose
    }

    day_folders = ParseList(pars)
    WriteCondorFiles(day_folders)
    WriteCondorScripts(day_folders, pars)
    SubmitJobs(day_folders)

if __name__ == '__main__':
    main()