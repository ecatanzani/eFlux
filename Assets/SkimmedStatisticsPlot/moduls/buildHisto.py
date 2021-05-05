import os
import argparse
import numpy as np
from tqdm import tqdm
from ROOT import TFile, TTree
import matplotlib.pyplot as plt
from datetime import date, timedelta

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

def fillStats(files: list, pars: dict, opts: argparse.Namespace) -> dict:

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

    for file in tqdm(files):
        rootfilename = pars['farmAddress'] + "/" + file
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
            status = True
            for erange in tmpstats:
                if tmpstats[erange] == -1:
                    status = False
                    break
            if status:
                evts_20_100.append(tmpstats['20_100'])
                evts_100_250.append(tmpstats['100_250'])
                evts_250_500.append(tmpstats['250_500'])
                evts_500_1.append(tmpstats['500_1'])
                evts_1_5.append(tmpstats['1_5'])
                evts_5.append(tmpstats['g5'])
                totevents = evts_20_100[-1]+evts_100_250[-1]+evts_250_500[-1]+evts_500_1[-1]+evts_1_5[-1]+evts_5[-1]
                evts_20_100.append(totevents)
            else:
                evts_20_100.append(0)
                evts_100_250.append(0)
                evts_250_500.append(0)
                evts_500_1.append(0)
                evts_1_5.append(0)
                evts_5.append(0)
                evts.append(0)

            for erange in tmpstats:
                tmpstats[erange] = -1
            timebins.append(sdate)
            sdate = tmpdate
    
    return {'timebins': timebins, 'evts': evts, 'evts_20_100': evts_20_100, 'evts_100_250': evts_100_250, 'evts_250_500': evts_250_500, 'evts_500_1': evts_500_1, 'evts_1_5': evts_1_5, 'evts_5': evts_5}

def buildHisto(opts: argparse.Namespace, files: list, pars: dict):
    
    stats = fillStats(files, pars, opts)
    
    # Plot
    plt.style.use('seaborn-whitegrid')
    plt.plot(stats['timebins'], stats['evts'], marker='p', label="all energy", color="dimgray", linestyle='dashed', linewidth=2, markersize=12)
    plt.plot(stats['timebins'], stats['evts_20_100'], marker='o', label="20 - 100 GeV", color="cornflowerblue", linestyle='dashed', linewidth=2, markersize=12)
    plt.plot(stats['timebins'], stats['evts_100_250'], marker='^', label="100 - 250 GeV", color="darkorange", linestyle='dashed', linewidth=2, markersize=12)
    plt.plot(stats['timebins'], stats['evts_250_500'], marker='*', label="250 - 500 GeV", color="forestgreen", linestyle='dashed', linewidth=2, markersize=12)
    plt.plot(stats['timebins'], stats['evts_500_1'], marker='|', label="0.5 - 1 TeV", color="crimson", linestyle='dashed', linewidth=2, markersize=12)
    plt.plot(stats['timebins'], stats['evts_1_5'], marker='d', label="1 - 5 TeV", color="blueviolet", linestyle='dashed', linewidth=2, markersize=12)
    plt.plot(stats['timebins'], stats['evts_5'], marker='s', label="> 5 TeV", color="saddlebrown", linestyle='dashed', linewidth=2, markersize=12)

    #plt.legend(numpoints=1)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., numpoints=1)
    plt.yscale('log')
    plt.ylim(10, 1e+6)
    plt.ylabel("counts")
    plt.savefig(opts.output)