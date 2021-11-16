import streamlit as st
from pathlib import Path
from datetime import datetime
import subprocess


def getDataFiles(start_date: datetime, end_date: datetime, farm_address: str, farm_skimmed_path: str) -> list:

    st.info("**Searching skimmed data files on XROOTD...**")
    st.info("This process may require some minutes accordingly to the selected time window ...")

    # Crate output data file list
    dList = []
    start = False

    # Get stage 0 dirs --> /FM/skim/6.0.0/v2/
    getDataDirsCommand = f"xrdfs {farm_address} ls {farm_skimmed_path}"
    dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
    dataDirs = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

    # Get stage 1 dirs --> /FM/skim/6.0.0/v2/YearOfAcquisition/
    for dir_st1 in dataDirs:
        if "20" in dir_st1 and "statistics" not in dir_st1:

            # Extract year
            year = int(dir_st1[dir_st1.rfind('/')+1:])
            if year < start_date.year or year > end_date.year:
                continue

            getDataDirsCommand = f"xrdfs {farm_address} ls {dir_st1}"
            dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
            dataDirs_st1 = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

            # Get stage 2 dirs --> /FM/skim/6.0.0/v2/YearOfAcquisition/MonthOfAcquisition
            for dir_st2 in dataDirs_st1:

                month = int(dir_st2[dir_st2.rfind('/')+1:])
                if not start and month < start_date.month:
                    continue
                elif start and month > end_date.month:
                    break

                getDataDirsCommand = f"xrdfs {farm_address} ls {dir_st2}"
                dataDirsOut = subprocess.run(getDataDirsCommand, shell=True, check=True, stdout=subprocess.PIPE)
                dataDirs_st2 = str.split(dataDirsOut.stdout.decode('utf-8').rstrip(), '\n')

                # Get ROOT data file
                for data_elm in dataDirs_st2:
                    if data_elm.endswith('.root') and "data_photon" not in data_elm:
                        day = int(data_elm[data_elm.rfind('/')+1:data_elm.find('_')])
                        if not start and day < start_date.day:
                            continue
                        elif not start and day >= start_date.day and day < end_date.day:
                            start = True
                        elif start and day > end_date.day:
                            break
                        
                        print(f"Adding file to the list: {data_elm}")
                        dList.append(data_elm)

    return dList
