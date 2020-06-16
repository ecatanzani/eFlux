import subprocess
from createFinalLists import createOutDataFile
import sys

def createDATAlist(pars, opts):
    dList = []

    # Get stage 0 dirs --> /FM/FlightData/2A/
    getDataDirsCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], pars['data_XRDFS_path'])
    if opts.verbose:
        print('Executing XRDFS command: {}'.format(getDataDirsCommand))
    dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
    dataDirs = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')
    
    # Get stage 1 dirs --> /FM/FlightData/2A/DayOfAcquisition/
    for dir_st1 in dataDirs:
        # Select interesting data dirs, that starts with the run year (2015, 2016, ...)
        if "2A/20" in dir_st1:
            getDataDirsCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], dir_st1)
            if opts.verbose:
                print('Executing XRDFS command: {}'.format(getDataDirsCommand))
            dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
            dataDirs_st1 = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

            # Get stage 2 dirs --> /FM/FlightData/2A/DayOfAcquisition/DataDir
            for dir_st2 in dataDirs_st1:
                getDataDirsCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], dir_st2)
                if opts.verbose:
                    print('Executing XRDFS command: {}'.format(getDataDirsCommand))
                dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
                dataDirs_st2 = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

                # Get ROOT data file
                for data_elm in dataDirs_st2:
                    if data_elm.endswith('.root'):
                        dataTmpCompletePath = pars['farmAddress'] + data_elm
                        dList.append(dataTmpCompletePath)
    
    with open(createOutDataFile(opts, pars['data_eMin'], pars['data_eMax']), "w") as outList:
        for file in dList:
            outPath = pars['farmAddress'] + file + "\n"
            outList.write(outPath)
