import streamlit as st
import pandas as pd
import altair as alt

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

def build_histos(rdf: pd.DataFrame):
    
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
