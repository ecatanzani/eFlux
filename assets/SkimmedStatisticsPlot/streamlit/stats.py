from stmoduls import app_settings
from search import getDataFiles 
from rdf import buildRDF, readRDF
from histos import build_histos
import streamlit as st
import pandas as pd

def main():
    
    start_date, end_date, source = app_settings()

    if source['data_storage_opt'] == 'XROOTD':
        start = st.sidebar.button('Search on XROOTD')
        if start:
            files = getDataFiles(start_date, end_date, source['farm_address'], source['farm_skimmed_path'])
            rdf = buildRDF(files, source['farm_address']) 
            build_histos(rdf)   
    else:
        if source['local_dir']:
            start = st.sidebar.button('Search local DataSet')
            if start:
                rdf = readRDF(source['local_dir'], start_date, end_date)
                build_histos(rdf)


if __name__ == "__main__":
    main()
