import subprocess
import sys
import os


def createDATAlist(pars, opts, sYear = 2016, eYear = 2020):
    if sYear > eYear:
        print("ERROR: Start year could not be bigger respect to end year")
        sys.exit()

    # Create data dictionary
    dataDict = {str(year) + "_" + str(month) : [] for year in range(sYear, eYear+1) for month in range(1, 13)}
    # Initialize data file counter
    counter = 0

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
            # Extract year
            s_year_idx = dir_st1.find("2A/20") + 3
            year = int(dir_st1[s_year_idx : s_year_idx + 4])
            # Extract month
            month = int(dir_st1[s_year_idx + 4 : s_year_idx + 6])

            if year < sYear or year > eYear:
                continue

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
                        dataDict[str(year) + "_" + str(month)].append(data_elm)
                        counter += 1
    
    if opts.verbose:
        print('{} data files have been read...'.format(counter))

    if opts.collapse:
        with open("dataFileList.txt", "w") as outList:
            for key in dataDict:
                for data_file in dataDict[key]:
                    outPath = pars['farmAddress'] + data_file + "\n"
                    outList.write(outPath)
    else:
        listPath = "dataList"
        os.mkdir(listPath)
        for year in [i for i in range(sYear, eYear +1)]:
            for month in range (1, 13):
                listPath = "dataList" + "/" + str(year) + "_" + str(month)
                os.mkdir(listPath)
                listPath += "/dataFileList.txt"
                with open(listPath, "w") as outList:
                    for data_file in dataDict[str(year) + "_" + str(month)]:
                        outPath = pars['farmAddress'] + data_file + "\n"
                        outList.write(outPath)
