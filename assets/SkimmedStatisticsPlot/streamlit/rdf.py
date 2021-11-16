import streamlit as st
import pandas as pd
from ROOT import TFile
from datetime import datetime, date, timedelta

def testROOTFile(path: str) -> bool:
    _tmp_file = TFile.Open(path)
    if _tmp_file and not _tmp_file.IsOpen():
        return False
    elif _tmp_file and _tmp_file.IsOpen() and _tmp_file.IsZombie():
        _tmp_file.Close()
        return False
    elif _tmp_file and _tmp_file.IsOpen() and _tmp_file.TestBit(TFile.kRecovered):
        _tmp_file.Close()
        return False
    else:
        _tmp_file.Close()
        return True

def getdate(file: str) -> date:
    year = int(file[file.rfind('/')-7:file.rfind('/')-3])
    month = int(file[file.rfind('/')-2:file.rfind('/')])
    day = int(file[file.rfind('/')+1:file.rfind('/')+3])
    return date(year, month, day)

def getdatelist(dstart: date, dend: date) -> list:
    delta = dend - dstart
    dlist = []
    for day in range(delta.days+1):
        dlist.append(dstart + timedelta(days=day))
    return dlist

def getStatsFromFile(rootfilename: str) -> int:
    tmpfile = TFile.Open(rootfilename, "READ")
    mytree = tmpfile.Get("CollectionTree")
    stats = mytree.GetEntries()
    tmpfile.Close()
    return stats

def fillStats(files: list, farm_address: str) -> dict:

    files.sort()
    # Get start & end dates
    sdate = getdate(files[0])

    # Build events arrays
    evts = []
    evts_20_100 = []
    evts_100_250 = []
    evts_250_500 = []
    evts_500_1 = []
    evts_1_5 = []
    evts_5 = []
    timebins = []

    tmpstats = {'20_100': -1, '100_250': -1, '250_500': -1, '500_1': -1, '1_5': -1, 'g5': -1}
    
    perc_complete = 0.
    step = 1./len(files)
    bar = st.progress(perc_complete)

    for file in files:
        
        rootfilename = f"{farm_address}/{file}"
        print(f"Reading file ... [{rootfilename}]")
        tmpdate = getdate(file)
        if tmpdate==sdate:
            if testROOTFile(rootfilename):
                filestats = getStatsFromFile(rootfilename)
                if "002_010" in file:
                    tmpstats['20_100'] = filestats
                elif "010_025" in file:
                    tmpstats['100_250'] = filestats
                elif "025_050" in file:
                    tmpstats['250_500'] = filestats
                elif "050_100" in file:
                    tmpstats['500_1'] = filestats
                elif "100_500" in file:
                    tmpstats['1_5'] = filestats
                elif "500_000" in file:
                    tmpstats['g5'] = filestats
                
        else:

            for erange in tmpstats:
                if tmpstats[erange] != -1:
                    if "20_100" in erange:
                        evts_20_100.append(tmpstats['20_100'])
                    elif "100_250" in erange:
                        evts_100_250.append(tmpstats['100_250'])
                    elif "250_500" in erange:
                        evts_250_500.append(tmpstats['250_500'])
                    elif "500_1" in erange:
                        evts_500_1.append(tmpstats['500_1'])
                    elif "1_5" in erange:
                        evts_1_5.append(tmpstats['1_5'])
                    elif "g5" in erange:
                        evts_5.append(tmpstats['g5'])

                else:
                    if "20_100" in erange:
                        evts_20_100.append(0)
                    elif "100_250" in erange:
                        evts_100_250.append(0)
                    elif "250_500" in erange:
                        evts_250_500.append(0)
                    elif "500_1" in erange:
                        evts_500_1.append(0)
                    elif "1_5" in erange:
                        evts_1_5.append(0)
                    elif "g5" in erange:
                        evts_5.append(0)

            totevents = evts_20_100[-1]+evts_100_250[-1]+evts_250_500[-1]+evts_500_1[-1]+evts_1_5[-1]+evts_5[-1]
            evts.append(totevents)


            for erange in tmpstats:
                tmpstats[erange] = -1
            timebins.append(sdate)
            sdate = tmpdate
            
            if testROOTFile(rootfilename):
                filestats = getStatsFromFile(rootfilename)
                if "002_010" in file:
                    tmpstats['20_100'] = filestats
                elif "010_025" in file:
                    tmpstats['100_250'] = filestats
                elif "025_050" in file:
                    tmpstats['250_500'] = filestats
                elif "050_100" in file:
                    tmpstats['500_1'] = filestats
                elif "100_500" in file:
                    tmpstats['1_5'] = filestats
                elif "500_000" in file:
                    tmpstats['g5'] = filestats
        
        # Update progress bar         
        perc_complete += step
        bar.progress(round(perc_complete, 1))

    return {'date': timebins, 'evts': evts, 'evts_20_100': evts_20_100, 'evts_100_250': evts_100_250, 'evts_250_500': evts_250_500, 'evts_500_1': evts_500_1, 'evts_1_5': evts_1_5, 'evts_5': evts_5}

def buildRDF(files: list, farm_address: str) -> pd.DataFrame:
    rdf = pd.DataFrame(fillStats(files, farm_address))
    rdf.to_csv('skimmed_data_files_rdf.csv', index=False)
    return rdf

def readRDF(path: str, start_date: datetime, end_date: datetime) -> pd.DataFrame:
    rdf = pd.read_csv(path)
    rdf['date'] = pd.to_datetime(rdf['date']).dt.date
    mask = (rdf['date'] >= start_date) & (rdf['date'] <= end_date)
    return rdf.loc[mask]



