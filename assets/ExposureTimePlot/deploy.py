from argparse import ArgumentParser
from tqdm import tqdm
from datetime import date
import subprocess
import os


def ExtractDateFromFile(file):
    year = int(file[file.rfind('_OBS_')+5:file.rfind('_OBS_')+9])
    month = int(file[file.rfind('_OBS_')+9:file.rfind('_OBS_')+11])
    day = int(file[file.rfind('_OBS_')+11:file.rfind('_OBS_')+13])
    return date(year, month, day)

def ParseInputFiles(input_files):
    with open(input_files, 'r') as inputlist:
        filelist = inputlist.read().splitlines()
    date = None
    day_file_list = []
    day_info = {}
    for file in tqdm(filelist, desc='Parsing input list...'):
        tmp_date =  ExtractDateFromFile(file)
        if date is None:
            date = tmp_date
        else:
            if date != tmp_date:
                day_info[date] = day_file_list
                date = tmp_date
                day_file_list = []
        day_file_list.append(file)
    if day_file_list:
        day_info[date] = day_file_list
    return day_info 

def SplitDaysInfo(days_info, days):
    days_info_split = []
    tmp_info_integrate = {}
    for idx, tmp_date in enumerate(days_info):
        tmp_info_integrate[tmp_date] = days_info[tmp_date]
        if not (idx+1)%days:
            days_info_split.append(tmp_info_integrate)
            tmp_info_integrate = {}

    if len(tmp_info_integrate):
        days_info_split.append(tmp_info_integrate)
    return days_info_split

def create_file_list(days_info, job_folder, verbose):
    dataListPath = '{}/dataList.txt'.format(job_folder)

    try:
        with open(dataListPath, 'w') as datalist:
            for tmp_date in days_info:
                for file in days_info[tmp_date]:
                    datalist.write('{}\n'.format(file))
    except OSError:
        print('ERROR creating data file list in: {}'.format(job_folder))
        raise
    else:
        if verbose:
            print('Data file list created in: {}'.format(job_folder))

def create_condor_submit_file(job_folder, verbose):
    outputPath = '{}/output.log'.format(job_folder)
    logPath = '{}/output.clog'.format(job_folder)
    errPath = '{}/output.err'.format(job_folder)
    bashScriptPath = '{}/script.sh'.format(job_folder)
    subFilePath = '{}/cndr.sub'.format(job_folder)

    try:
        with open(subFilePath, 'w') as outSub:
            outSub.write("universe = vanilla\n")
            outSub.write('executable = {}\n'.format(bashScriptPath))
            outSub.write('output = {}\n'.format(outputPath))
            outSub.write('error = {}\n'.format(errPath))
            outSub.write('log = {}\n'.format(logPath))
            outSub.write("ShouldTransferFiles = YES\n")
            outSub.write("WhenToTransferOutput = ON_EXIT\n")
            outSub.write("queue 1")
    except OSError:
        print('ERROR creating HTCondor sub file in: {}'.format(job_folder))
        raise
    else:
        if verbose:
            print('HTCondor sub file created in: {}'.format(job_folder))

def create_bash_script(job_folder, python_script, verbose):
    bashScriptPath = '{}/script.sh'.format(job_folder)
    tmpOutDir = '{}/outFiles'.format(job_folder)
    dataListPath = '{}/dataList.txt'.format(job_folder)
    
    try:
        with open(bashScriptPath, "w") as outScript:
            
            outScript.write("#!/usr/bin/env bash\n")
            outScript.write("source /cvmfs/dampe.cern.ch/centos7/etc/setup.sh\n")
            outScript.write("dampe_init trunk\n")
            outScript.write('mkdir {}\n'.format(tmpOutDir))
            outScript.write('python {} -i {} -o {} -v\n'.format(python_script, dataListPath, os.path.join(tmpOutDir, "out")))

    except OSError:
        print('ERROR creating HTCondor bash script file in: {}'.format(job_folder))
        raise
    else:
        if verbose:
            print('HTCondor bash script file created in: {}'.format(job_folder))

def create_job_folder(days_info, python_script, condor_output, idx, verbose):
    job_folder = os.path.join(condor_output, 'job_'+str(idx))
    if not os.path.exists(job_folder):
        os.makedirs(job_folder)
    
    # Create data file list
    create_file_list(days_info, job_folder, verbose)
    # Create HTCondor submit file
    create_condor_submit_file(job_folder, verbose)
    # Create bash script file
    create_bash_script(job_folder, python_script, verbose)
    
def submit_job(condor_output, idx, verbose):
    job_folder = os.path.join(condor_output, 'job_'+str(idx))
    subFilePath = '{}/cndr.sub'.format(job_folder)
    subprocess.run(['condor_submit -name sn-02.cr.cnaf.infn.it -spool {}'.format(subFilePath)], shell=True, check=True)

def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Flight DATA statistics - Condor Engine")
    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='input DATA list')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output HTCondor directory')
    parser.add_argument("-d", "--days", type=int, dest='days',
                        const=10, nargs='?', help='Days to process in a single job')
    parser.add_argument("-s", "--script", type=str,
                        dest='script', help='Analysis python script')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)

    # Parse input file list
    days_info = ParseInputFiles(opts.input)
    # Split dictionary into days for jobs submission
    splitted_days_info = SplitDaysInfo(days_info, opts.days)

    # Prepare job folders
    for idx, days_info in enumerate(splitted_days_info):
        create_job_folder(days_info, opts.script, opts.output, idx, opts.verbose)
        submit_job(opts.output, idx, opts.verbose)

if __name__ == '__main__':
    main()