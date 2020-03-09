#!/usr/bin/env python3
from datetime import datetime
from argparse import ArgumentParser
import json
import subprocess

start = datetime.now()

def parseConfigFile(opts):
    dConfig = {'farmAddress': "", 'simu_XRDFS_path': "", 'geometry': "",
               'simu_eMin': 0, 'simu_eMax': 0, 'particle': "", 'jSet': ""}
    
    config_params = []
    custom_farm = False
    custom_set = False

    if opts.farm:
        farmAddress = opts.farm
        dConfig['farmAddress'] = farmAddress
        custom_farm = True
    if opts.input:
        jSet = opts.input
        dConfig['jSet'] = jSet
        custom_set = True
    
    
    with open("crawlerConfig.txt", "r") as _config:
        for line in _config:
            for word in line.split():
                config_params.append(word)

    for idx, word in enumerate(config_params):
        if word == "farmAddress":
            if not custom_farm:
                dConfig['farmAddress'] = config_params[idx+1]
        if word == "simu_XRDFS_path":
            dConfig['simu_XRDFS_path'] = config_params[idx+1]
        if word == "geometry":
            dConfig['geometry'] = config_params[idx+1]
        if word == "simu_eMin":
            dConfig['simu_eMin'] = float(config_params[idx+1])
        if word == "simu_eMax":
            dConfig['simu_eMax'] = float(config_params[idx+1])
        if word == "particle":
            dConfig['particle'] = config_params[idx+1]
        if word == "jSet":
            if not custom_set:
                dConfig['jSet'] = config_params[idx+1]
    
    return dConfig

def createOutFile(opts, simu_eMin, simu_eMax, particle):
    if not opts.input:
        finalPath = "simuFileList_"
        finalPath += particle
        finalPath += "_"
        finalPath += str(simu_eMin)
        finalPath += "_"
        finalPath += str(simu_eMax)
        finalPath += ".txt"
        return finalPath
    else:
        return opts.output


def main(args=None):

    parser = ArgumentParser(
        usage="Usage: %(prog)s [options]", description="reco MC files explorer")

    parser.add_argument("-i", "--input", type=str,
                        dest='input', help='json data storage')
    parser.add_argument("-f", "--farm", type=str,
                        dest='farm', help='MC ROOT input farm')
    parser.add_argument("-o", "--output", type=str,
                        dest='output', help='output file list')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)
    
    # Get dictionary from config file parsing
    pars = parseConfigFile(opts)

    # Parse data json
    with open(pars['jSet'], 'rb') as set_:
        dataSets = json.load(set_)

        # Create out file list
        dList = []
        for elm in dataSets:
            if elm["geometry"] == pars['geometry']:
                if elm["particle"] == pars['particle']:
                    if int(elm["eMin"]) >= pars['simu_eMin']:
                        if int(elm["eMax"]) <= pars['simu_eMax']:
                            tmpSetPath = pars['simu_XRDFS_path'] + "v" + \
                                pars['geometry'] + "/" + elm["dSet"]
                            dList.append(tmpSetPath)

    with open(createOutFile(opts, pars['simu_eMin'], pars['simu_eMax'], pars['particle']), "w") as outList:
        for set_ in dList:
            tmpCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], set_)
            if opts.verbose:
                print(tmpCommand)
            tmpOut = subprocess.run(
                tmpCommand, shell=True, check=True, stdout=subprocess.PIPE)
            setList = str.split(tmpOut.stdout.decode('utf-8').rstrip(), '\n')
            for file in setList:
                outPath = pars['farmAddress'] + file + "\n"
                outList.write(outPath)


if __name__ == "__main__":
    main()
    print(datetime.now()-start)
