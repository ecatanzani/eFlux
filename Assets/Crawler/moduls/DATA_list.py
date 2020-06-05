import subprocess
from createFinalLists import createOutDataFile


def createDATAlist(pars, opts):
    dList = []
    if pars['geometry'] == "6r0p0":
        geoPars = "6.0.0"
    if pars['geometry'] == "5r4p0":
        geoPars = "5.4.0"
    dataPath = pars['data_XRDFS_path'] + geoPars

    getDataDirsCommand = 'xrdfs {} ls {}'.format(
        pars['farmAddress'], dataPath)
    if opts.verbose:
        print('Executing XRDFS command: {}'.format(getDataDirsCommand))
    dataDirsOut = subprocess.run(
        getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
    dataDirs = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')
    for elm in dataDirs:
        tmpDataCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], elm)
        tmpData = subprocess.run(
            tmpDataCommand, shell=True, check=True, stdout=subprocess.PIPE)
        dataInDir = str.split(
            tmpData.stdout.decode('utf-8').rstrip(), '\n')
        for sData in dataInDir:
            if sData.endswith('.root'):
                dataTmpCompletePath = pars['farmAddress'] + sData
                dList.append(dataTmpCompletePath)
    with open(createOutDataFile(opts, pars['data_eMin'], pars['data_eMax']), "w") as outList:
        for file in dList:
            outPath = pars['farmAddress'] + file + "\n"
            outList.write(outPath)
