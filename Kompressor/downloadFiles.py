import subprocess
from argparse import ArgumentParser


def main(args=None):
    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="DATA/MC file downloader")
    parser.add_argument("-l", "--list", type=str,
                        dest='list', help='input file list')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output directory')
    parser.add_argument("-s", "--scp", type=str,
                        dest='scp', help='ssh config (if scp protol has been choosen)')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)

    with open(opts.list, 'r') as inputlist:
        filelist = inputlist.read().splitlines()
    
    for idx, file in enumerate(filelist):
        if opts.verbose:
            print(f"Downloading file: {file}")
        if opts.scp:
            file_name = f"{file[file.rfind('/')+1:file.rfind('.root')]}_{idx}.root"
            subprocess.run(f"scp {opts.scp}:{file} {opts.output}/{file_name}", shell=True, check=True, stdout=subprocess.PIPE)
        else:
            subprocess.run(f"xrdcp {file} {opts.output}", shell=True, check=True, stdout=subprocess.PIPE)

if __name__ == '__main__':
    main()
