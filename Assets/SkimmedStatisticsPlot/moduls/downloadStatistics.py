from pathlib import Path
import subprocess
import argparse
import sys
import os


def getStatisticFiles(config: dict, opts: argparse.Namespace) -> list:

    # Crate output data file list
    dList = []
    years = []
    counters = []

    # Get stage 0 dirs --> /FM/skim/6.0.0/v2/
    getDataDirsCommand = f"xrdfs {config['farmAddress']} ls {config['data_XRDFS_skimmed_path']}"
    if opts.verbose:
        print(f"Executing XRDFS command: {format(getDataDirsCommand)}")
    dataDirsOut = subprocess.run(
        getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
    dataDirs = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

    # Get stage 1 dirs --> /FM/FlightData/2A/YearOfAcquisition/
    for dir_st1 in dataDirs:
        if "20" in dir_st1 and "statistics" not in dir_st1:

            # Extract year
            year = int(dir_st1[dir_st1.rfind('/')+1:])
            if year not in years:
                years.append(year)
                year_data_idx = len(years)-1
                counters.append(0)

            getDataDirsCommand = f"xrdfs {config['farmAddress']} ls {dir_st1}"
            if opts.verbose:
                print(f"Executing XRDFS command: {getDataDirsCommand}")
            dataDirsOut = subprocess.run(
                getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
            dataDirs_st1 = str.split(
                dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

            # Get stage 2 dirs --> /FM/FlightData/2A/YearOfAcquisition/MonthOfAcquisition
            for dir_st2 in dataDirs_st1:
                getDataDirsCommand = f"xrdfs {config['farmAddress']} ls {dir_st2}"
                if opts.verbose:
                    print(f"Executing XRDFS command: {getDataDirsCommand}")
                dataDirsOut = subprocess.run(
                    getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
                dataDirs_st2 = str.split(
                    dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

                # Get ROOT data file
                for data_elm in dataDirs_st2:
                    if data_elm.endswith('.root') and "data_photon" not in data_elm:
                        dList.append(data_elm)
                        counters[year_data_idx] += 1

    if opts.verbose:
        print(f"{sum(counters)} data files have been read...")
        for year_idx, year in enumerate(years):
            print(f"{counters[year_idx]} data files found in {year} folder")

    return dList
