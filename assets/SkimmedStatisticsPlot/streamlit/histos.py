import streamlit as st
import pandas as pd
import altair as alt
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.dates as mdates


from ROOT import TFile, TH1D
import sys

def single_histo(rdf: pd.DataFrame, erange: str, mycolor: str, meanbarcolor: str = 'red'):
    
    base =  alt.Chart(rdf)
    
    histo = base.mark_bar().mark_area(
        opacity=0.3,
        interpolate='step'
    ).encode(
        x=alt.X(f'{erange}:Q', bin=True),
        y=alt.Y('count()'),
        color=alt.value(mycolor)
    ).properties(title='Event Statistics').interactive()


    rule = base.mark_rule(color=meanbarcolor).encode(
        x=f'mean({erange}):Q',
        size=alt.value(3)
    )

    st.write(histo+rule)

def build_ROOT_histos(rdf: pd.DataFrame):

    out_file = TFile('histos.root', 'RECREATE')
    if not out_file.IsOpen():
        print('Could not write ROOT output file')
        sys.exit()
    
    h_stat = TH1D('h_stat', 'h_stat', 100, 150000, 300000)
    h_stat_20_100 = TH1D("h_stat_20_100", "data statistics - 20 GeV - 100 GeV", 100, 150000, 300000)
    h_stat_100_250 = TH1D("h_stat_100_250", "data statistics - 100 GeV - 250 GeV", 100, 10000, 20000)
    h_stat_250_500 = TH1D("h_stat_250_500", "data statistics - 250 GeV - 500 GeV", 100, 2000, 4000)
    h_stat_500_1 = TH1D("h_stat_500_1", "data statistics - 500 GeV - 1 TeV", 100, 100, 2000)
    h_stat_1_5 = TH1D("h_stat_1_5", "data statistics - 1 TeV - 5 TeV", 100, 0, 1000)
    h_stat_5 = TH1D("h_stat_5", "data statistics - 5 TeV", 50, 0, 100)

    for day_idx in range(len(rdf.index)):
        h_stat.Fill((rdf['evts'].values)[day_idx])
        h_stat_20_100.Fill((rdf['evts_20_100'].values)[day_idx])
        h_stat_100_250.Fill((rdf['evts_100_250'].values)[day_idx])
        h_stat_250_500.Fill((rdf['evts_250_500'].values)[day_idx])
        h_stat_500_1.Fill((rdf['evts_500_1'].values)[day_idx])
        h_stat_1_5.Fill((rdf['evts_1_5'].values)[day_idx])
        h_stat_5.Fill((rdf['evts_5'].values)[day_idx])

    h_stat.Write()
    h_stat_20_100.Write()
    h_stat_100_250.Write()
    h_stat_250_500.Write()
    h_stat_500_1.Write()
    h_stat_1_5.Write()
    h_stat_5.Write()

    out_file.Close()

def write_final_graph(rdf: pd.DataFrame):

    rcParams.update({'figure.autolayout': True})
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=3))

    plt.plot(rdf['date'], rdf['evts'], label="all energies", color="dimgray")
    plt.plot(rdf['date'], rdf['evts_20_100'], label="20 - 100 GeV", color="cornflowerblue")
    plt.plot(rdf['date'], rdf['evts_100_250'], label="100 - 250 GeV", color="darkorange")
    plt.plot(rdf['date'], rdf['evts_250_500'], label="250 - 500 GeV", color="forestgreen")
    plt.plot(rdf['date'], rdf['evts_500_1'], label="0.5 - 1 TeV", color="crimson")
    plt.plot(rdf['date'], rdf['evts_1_5'], label="1 - 5 TeV", color="blueviolet")
    plt.plot(rdf['date'], rdf['evts_5'], label="> 5 TeV", color="saddlebrown")
    
    plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
    plt.yscale('log')
    plt.ylim(10, 1e+6)
    plt.ylabel("counts/day")

    plt.gcf().autofmt_xdate()
    plt.savefig('stats.pdf')

def build_histos(rdf: pd.DataFrame):

    rdf["date"] = pd.to_datetime(rdf["date"])
    
    gray = '#7f7f7f'
    blue ='#1f77b4'
    orange = '#ff7f0e'
    green = '#2ca02c'
    red = '#d62728'
    brown = '#8c564b'
    violet = '#9467bd'

    color_range = [gray, blue, orange, green, red, violet, brown]

    all_line = alt.Chart(rdf).transform_fold(
        ['evts', 'evts_20_100', 'evts_100_250', 'evts_250_500', 'evts_500_1', 'evts_1_5', 'evts_5'],
    ).mark_line().encode(
        x='date:T',
        y=alt.Y("value:Q", scale=alt.Scale(type='log'), title='counts/day'),
        color=alt.Color('key:N', scale=alt.Scale(domain=['evts', 'evts_20_100', 'evts_100_250', 'evts_250_500', 'evts_500_1', 'evts_1_5', 'evts_5'], range=color_range, type="ordinal"))
    ).properties(title='Skimmed Data Statistics - Time Distribution').interactive()

    all_points = alt.Chart(rdf).transform_fold(
        ['evts', 'evts_20_100', 'evts_100_250', 'evts_250_500', 'evts_500_1', 'evts_1_5', 'evts_5'],
    ).mark_circle(
        size=80,
        filled=True,
        opacity=0.5
    ).encode(
        x='date:T',
        y=alt.Y("value:Q", scale=alt.Scale(type='log'), title='counts/day'),
        color=alt.Color('key:N', scale=alt.Scale(domain=['evts', 'evts_20_100', 'evts_100_250', 'evts_250_500', 'evts_500_1', 'evts_1_5', 'evts_5'], range=color_range, type="ordinal"))
    ).properties(title='Skimmed Data Statistics - Time Distribution').interactive()
    
    chart, histo = st.columns(2)
    with chart:
        st.write(all_line+all_points)
    with histo:
        single_histo(rdf, 'evts', gray)

    h_20_100, h_100_250, h_250_500 = st.columns(3)
    with h_20_100:
        single_histo(rdf, 'evts_20_100', blue)
    with h_100_250:
        single_histo(rdf, 'evts_100_250', orange)
    with h_250_500:
        single_histo(rdf, 'evts_250_500', green)

    h_500_1, h_1_5, h_5 = st.columns(3)
    with h_500_1:
        single_histo(rdf, 'evts_500_1', red, meanbarcolor='blue')
    with h_1_5:
        single_histo(rdf, 'evts_1_5', violet)
    with h_5:
        single_histo(rdf, 'evts_5', brown)

    build_ROOT_histos(rdf)
    write_final_graph(rdf)
