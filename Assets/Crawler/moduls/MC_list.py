import subprocess
import json
import os
import sys
from createFinalLists import createOutSimuFile


def createMClist(pars, opts):
    # Parse masks sets
    maskInputFile = os.getcwd()
    mask_sIdx = maskInputFile.find('/Crawler')
    maskInputFile = maskInputFile[0:mask_sIdx+8] + "/setMask.txt"
    try:
        dListMask = open(maskInputFile).read().splitlines()
    except OSError:
            print('ERROR reading input dataSet mask file in: {}'.format(maskInputFile))
            raise
    # Parse data json
    with open(pars['jSet'], 'rb') as set_:
        dataSets = json.load(set_)
        # Create out file list
        dList = []
        setPowerLowConfig = []
        for elm in dataSets:
            if elm["dSet"] in dListMask:
                continue
            if elm["geometry"] == pars['geometry']:
                if elm["particle"] == pars['particle']:
                    if int(elm["eMin"]) >= pars['simu_eMin']:
                        if int(elm["eMax"]) <= pars['simu_eMax']:
                            tmpSetPath = pars['simu_XRDFS_path'] + "v" + \
                                pars['geometry'] + "/" + elm["dSet"]
                            dList.append(tmpSetPath)
                            setPowerLowConfig.append(
                                [elm['dSet'], elm['eMin'], elm['eMax'], elm['pSpectrum']])

    finalListPath = createOutSimuFile(
        opts, pars['simu_eMin'], pars['simu_eMax'], pars['particle'])
    if opts.verbose:
        print('Writing output MC file list: {}'.format(finalListPath))
    with open(finalListPath, "w") as outList:
        for set_ in dList:
            tmpCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], set_)
            if opts.verbose:
                print('Executing XRDFS command: {}'.format(tmpCommand))
            tmpOut = subprocess.run(
                tmpCommand, shell=True, check=True, stdout=subprocess.PIPE)
            setList = str.split(tmpOut.stdout.decode('utf-8').rstrip(), '\n')
            for file in setList:
                outPath = pars['farmAddress'] + file + "\n"
                outList.write(outPath)

    currentWD = os.getcwd()
    assets_index = currentWD.find('/Assets')
    powerLowConfigPath = currentWD[0:assets_index] + \
        "/config/powerlowMCsets.txt"
    if opts.verbose:
        print('Writing output MC power law config file: {}'.format(
            powerLowConfigPath))
    with open(powerLowConfigPath, "w") as pLawConfig:
        pLawConfig.write(
            "######### Power Law DataSets Config File #########\n\n")
        pLawConfig.write("SetName\t eMin\t eMax\t generationPawerLawIndex\n")
        for _set in setPowerLowConfig:
            tmpString = "Set:"
            for elm in _set:
                tmpString += "\t" + str(elm)
            tmpString += "\n"
            pLawConfig.write(tmpString)
