#!/usr/bin/env python3
from datetime import datetime
from argparse import ArgumentParser
import json
import subprocess

start = datetime.now()

def parseConfigFile(opts):
    dConfig = {'farmAddress': "", 'simu_XRDFS_path': "", 'data_XRDFS_path': "", 'geometry': "", 'simu_eMin': 0, 'simu_eMax': 0, 'data_eMin': 0, 'data_eMax': 0, 'particle': "", 'jSet': ""}
    
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
        if word == "data_XRDFS_path":
            dConfig['data_XRDFS_path'] = config_params[idx+1]
        if word == "geometry":
            dConfig['geometry'] = config_params[idx+1]
        if word == "simu_eMin":
            dConfig['simu_eMin'] = float(config_params[idx+1])
        if word == "simu_eMax":
            dConfig['simu_eMax'] = float(config_params[idx+1])
        if word == "data_eMin":
            dConfig['data_eMin'] = float(config_params[idx+1])
        if word == "data_eMax":
            dConfig['data_eMax'] = float(config_params[idx+1])
        if word == "particle":
            dConfig['particle'] = config_params[idx+1]
        if word == "jSet":
            if not custom_set:
                dConfig['jSet'] = config_params[idx+1]
    
    return dConfig

def createOutSimuFile(opts, simu_eMin, simu_eMax, particle):
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

def createOutSimuFile(opts, data_eMin, data_eMax):
    if not opts.input:
        finalPath = "dataFileList_"
        finalPath += str(data_eMin)
        finalPath += "_"
        finalPath += str(data_eMax)
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
    parser.add_argument("-d", "--data", dest='data', default=False,
                        action='store_true', help='get DATA file list')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False,
                        action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)
    
    # Get dictionary from config file parsing
    pars = parseConfigFile(opts)

    if opts.data:
        dList=[]
        if pars['geometry'] == "6r0p0":
            geoPars = "6.0.0"
        if pars['geometry'] == "5r4p0":
            geoPars = "5.4.0"
        dataPath = pars['data_XRDFS_path'] + geoPars

        getDataDirsCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], dataPath)
        if opts.verbose:
            print(getDataDirsCommand)
        dataDirsOut = subprocess.run(
            getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
        dataDirs = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')
        for elm in dataDirs:
            tmpDataCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], elm)
            tmpData = subprocess.run(
                tmpDataCommand, shell=True, check=True, stdout=subprocess.PIPE)
            dataInDir = str.split(tmpData.stdout.decode('utf-8').rstrip(), '\n')
            for sData in dataInDir:
                if sData.endswith('.root'):
                    dataTmpCompletePath = pars['farmAddress'] + sData
                    dList.append(dataTmpCompletePath)
        with open(createOutDataFile(opts, pars['data_eMin'], pars['data_eMax']), "w") as outList:
            for file in dList:
                outPath = pars['farmAddress'] + file + "\n"
                outList.write(outPath)

    else:
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

        with open(createOutSimuFile(opts, pars['simu_eMin'], pars['simu_eMax'], pars['particle']), "w") as outList:
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
