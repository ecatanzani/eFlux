import os
import sys
import argparse
from datetime import date


def getdate(file: str) -> date:
    year = int(file[file.rfind('/')+1:file.rfind('/')+5])
    month = int(file[file.rfind('/')+6:file.rfind('/')+8])
    day = int(file[file.rfind('/')+9:file.rfind('/')+11])
    return date(year, month, day)

def getfilename(_date: date, working_dir: str = "stats"):
    year = str(_date.year)
    month = str("0" + str(_date.month)) if _date.month < 10 else str(_date.month)
    day = str("0" + str(_date.day)) if _date.day < 10 else str(_date.day)
    return [ working_dir + "/" + year + "_" + month + "_" + day + "_data_" + erange + ".root.stats" for erange in ["002_010", "010_025", "025_050", "050_100", "100_500", "500_000"]]

def recover_files(_date: date, working_dir: str = "stats"):
    filelist = getfilename(_date, working_dir)
    for file in filelist:
        if os.path.exists(file):
            os.remove(file)
        with open(file, 'w') as wfile:
            wfile.write(f"Number of skimmed events in {file[file.rfind('/')+1:]}: 0")

def checkFiles(opts: argparse.Namespace, working_dir: str = "stats"):
    if opts.verbose:
        print(f"Checking downloaded files in {working_dir}")
        
    files = [working_dir + "/" + file for file in os.listdir(working_dir) if file.endswith(".root.stats")]
    files.sort()
    
    if len(files)%6 != 0:
        print("ERROR downloading files: number of files mismatch")

    file_status = {'20_100': False, '100_250': False, '250_500': False, '500_1': False, '1_5': False, 'g5': False}
    sdate = getdate(files[0])

    for file in files:
        tmpdate = getdate(file)
        if tmpdate==sdate:
            if "002_010" in file:
                file_status['20_100'] = True
            elif "010_025" in file:
                file_status['100_250'] = True
            elif "025_050" in file:
                file_status['250_500'] = True
            elif "050_100" in file:
                file_status['500_1'] = True
            elif "100_500" in file:
                file_status['1_5'] = True
            elif "500_000" in file:
                file_status['g5'] = True
        else:
            
            status = True
            for erange in file_status:
                status *= file_status[erange]
            if status:
                for erange in file_status:
                    file_status[erange] = False
            else:
                print(f"ERROR number of files... recovering -> {file}")
                recover_files(sdate, working_dir)

            if "002_010" in file:
                file_status['20_100'] = True
            elif "010_025" in file:
                file_status['100_250'] = True
            elif "025_050" in file:
                file_status['250_500'] = True
            elif "050_100" in file:
                file_status['500_1'] = True
            elif "100_500" in file:
                file_status['1_5'] = True
            elif "500_000" in file:
                file_status['g5'] = True

            sdate = tmpdate