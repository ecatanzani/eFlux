import os
import datetime
from tqdm import tqdm
from argparse import ArgumentParser

def ExtractDateFromFile(file: str) -> datetime.date:
    year = int(file[file.rfind('/2A/')+4:file.rfind('/2A/')+8])
    month = int(file[file.rfind('/2A/')+8:file.rfind('/2A/')+10])
    day = int(file[file.rfind('/2A/')+10:file.rfind('/2A/')+12])
    return datetime.date(year, month, day)

def WiteListOnFile(list: list, list_date: datetime.date, output: str):

    month = str(list_date.month) if list_date.month > 9 else f"0{list_date.month}"
    day = str(list_date.day) if list_date.day > 9 else f"0{list_date.day}"
    newdir = f"{output}/{list_date.year}{month}{day}"
    os.mkdir(newdir)
    with open(f"{newdir}/dataDayList.txt", 'w') as outputlist:
        for file in list:
            outputlist.write(f"{file}\n")

def ParseList(list: str, output_dir: str, verbose: bool):

    with open(list, 'r') as inputlist:
        filelist = inputlist.read().splitlines()
    
    date = None
    day_file_list = []
    for file in tqdm(filelist, desc='Parsing input list...'):
        tmp_date =  ExtractDateFromFile(file)
        if date is None:
            date = tmp_date
        else:
            if date != tmp_date:
                WiteListOnFile(day_file_list, date, output_dir)
                date = tmp_date
                day_file_list = []
        day_file_list.append(file)
    
    if day_file_list:
        WiteListOnFile(day_file_list, date, output_dir)



def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="Flight DATA statistics")
    parser.add_argument("-l", "--list", type=str,
                        dest='list', help='input file list')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output directory')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)

    ParseList(opts.list, opts.output, opts.verbose)




if __name__ == '__main__':
    main()