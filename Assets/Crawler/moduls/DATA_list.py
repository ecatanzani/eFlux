import subprocess
import sys
import os


def createDATAlist(pars, opts):
    if pars['data_sYear'] > pars['data_eYear']:
        print("ERROR: Start year could not be bigger respect to end year")
        sys.exit()

    # Crate output data file list
    dList = []
    years = []
    counters = []

    # Get stage 0 dirs --> /FM/FlightData/2A/
    getDataDirsCommand = 'xrdfs {} ls {}'.format(pars['farmAddress'], pars['data_XRDFS_path'])
    if opts.verbose:
        print('Executing XRDFS command: {}'.format(getDataDirsCommand))
    dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
    dataDirs = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

    if opts.verbose:
        print('Collecting data from {} to {}'.format(pars['data_sYear'], pars['data_eYear']))

    # Get stage 1 dirs --> /FM/FlightData/2A/DayOfAcquisition/
    for dir_st1 in dataDirs:

        # Select interesting data dirs, that starts with the run year (2015, 2016, ...)
        if "2A/20" in dir_st1:
            # Extract year
            s_year_idx = dir_st1.find("2A/20") + 3
            year = int(dir_st1[s_year_idx : s_year_idx + 4])

            if year < pars['data_sYear'] or year > pars['data_eYear']:
                continue

            if year not in years:
                years.append(year)
                year_data_idx = len(years)-1
                counters.append(0)

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
                        dList.append(data_elm)
                        counters[year_data_idx] += 1
    
    if opts.verbose:
        print('{} data files have been read...'.format(sum(counters)))
        for year_idx, year in enumerate(years):
            print('{} data files found in {} folder'.format(counters[year_idx], year))

    if opts.output:
        data_list_path = opts.output
    else:
        data_list_path =  "dataFileList.txt"
    with open(data_list_path, "w") as outList:
        for elm in dList:
            outList.write(pars['farmAddress'] + elm + "\n")
