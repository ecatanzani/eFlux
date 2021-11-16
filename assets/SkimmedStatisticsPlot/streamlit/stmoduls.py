import streamlit as st
import datetime

def app_settings() -> tuple:
    version = 1.0
    train = "stable"
    st.set_page_config(layout="wide")
    st.write(f"""
    # DAMPE skimmed data statistics
    release v{version} - train {train}
    """)

    source = {
        "data_storage_opt": str(),
        "farm_address": str(),
        "farm_skimmed_path": str(),
        "local_dir": str()
    }

    st.sidebar.write('**Settings**')

    start_date = st.sidebar.date_input('Start date', value=datetime.datetime(2016, 1, 1), help='Select starting data')
    end_date = st.sidebar.date_input('End date', value=datetime.date.today(), help='Select ending data')

    st.sidebar.write('**DATA source**')
    source['data_storage_opt'] = st.sidebar.selectbox('Choose how to get skimmed data files', ("XROOTD", 'Use local dir'), index=1)
    if source['data_storage_opt'] == 'XROOTD':
        source['farm_address'] = st.sidebar.text_input('XROOTD DAMPE entrypoint:', "root://xrootd-dampe.cloud.ba.infn.it//")
        source['farm_skimmed_path'] = st.sidebar.text_input('XROOTD DAMPE skimmed data files:', "/FM/skim/6.0.0/v2/")
    else:
        source['local_dir'] = st.sidebar.file_uploader("Upload a DB", help='Pick a skimmed data DB', accept_multiple_files=False, type=['csv'])

    return (start_date, end_date, source)