#!/usr/bin/env python3
from datetime import datetime
from argparse import ArgumentParser
import json
import subprocess

start = datetime.now()

def createOutFile(opts,simu_eMin,simu_eMax,particle):
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
    
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="reco MC files explorer")

    parser.add_argument("-i", "--input", type=str, dest='input', help='json data storage')
    parser.add_argument("-f", "--farm", type=str, dest='farm', help='MC ROOT input farm')
    parser.add_argument("-o", "--output", type=str, dest='output', help='output file list')
    parser.add_argument("-v", "--verbose", dest='verbose', default=False, action='store_true', help='run in high verbosity mode')

    opts = parser.parse_args(args)

    # DataSets default params
    farmAddress = "root://xrootd-dampe.cloud.ba.infn.it//"
    simu_XRDFS_path = "/MC/reco/"
    geometry = "6r0p0"
    simu_eMin = 0
    simu_eMax = 1000
    particle = "e"
    jSet = "dataSets.json"

    ###############################
    
    if opts.farm:
        farmAddress = opts.farm
    
    if opts.input:
        jSet = opts.input
    
    # Parse data json
    with open(jSet,'rb') as set_:
        dataSets = json.load(set_)
        
        # Create out file list
        dList = []
        for elm in dataSets:
            if elm["geometry"]==geometry:
                if elm["particle"]==particle:
                    if int(elm["eMin"])>=simu_eMin:
                        if int(elm["eMax"])<=simu_eMax:
                            tmpSetPath = simu_XRDFS_path + "v" + geometry + "/" + elm["dSet"]
                            dList.append(tmpSetPath)

    with open(createOutFile(opts,simu_eMin,simu_eMax,particle),"w") as outList:                
        for set_ in dList:
            tmpCommand = 'xrdfs {} ls {}'.format(farmAddress,set_)
            if opts.verbose:
                print(tmpCommand)
            tmpOut = subprocess.run(tmpCommand, shell=True, check=True, stdout=subprocess.PIPE)
            setList = str.split(tmpOut.stdout.decode('utf-8').rstrip(),'\n')
            for file in setList:
                outPath = farmAddress + file + "\n"
                outList.write(outPath)

if __name__ == "__main__":
    main()
    print(datetime.now()-start)