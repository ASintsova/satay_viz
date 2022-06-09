from src.count_inserts import IntervalCounter
from src.compare_interval_counts import concatenate_results
import streamlit as st
import pandas as pd
from datetime import datetime as dt

@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')


def app():
    st.header("Counting transposons over custom intervals")
    st.markdown("#### Step 1: Upload the `bed` files you want to analyze")
    today = dt.now().strftime("%d-%m-%Y")
    input_files = st.file_uploader('Upload bed files with transposon coordinates', accept_multiple_files=True)
    st.markdown("#### Step 2: Choose the length of the interval and an offset. ")
    st.markdown("For example, if the interval length = 500, and offset = 50 "
                "The transposons will be first counted over intervals of 1:500, 501:1000..."
                "and then over intervals 50:550, 551:1050, ... etc.")
    c1, c2 = st.columns(2)
    interval_size = c1.number_input('Choose interval size', value=500, min_value=100)
    offset = c2.number_input("Choose offset", value=50, min_value=25)
    chr_sizes = 'chr_sizes.tsv'
    st.markdown("#### Step 3: Enter experiment name. This will be used for the name of final output file.")
    exp_name = st.text_input("Enter experiment name", value=f'{today}_{interval_size}int_{offset}offset')
    ic = IntervalCounter(input_files, chr_sizes, interval_size, offset, exp_name, '.', sep=' ', skiprows=1)
    st.markdown("#### Step 4: Count.")
    if st.button('Count'):
        if not input_files:
            st.markdown('No `bed` files provided')
            st.stop()
        st.markdown("Counting ...")
        ic.count()
        offset_counts = pd.concat([r.as_df().assign(offset=t.split('.')[0].split("_")[-1]) for t, r in ic.overlaps.items()])
        csv = convert_df(offset_counts)
        st.markdown("Counting complete")
        st.markdown("#### Step 5: Download and save your results. This file can be used with `Analysis` page")
        st.download_button(
            label="Download counts as CSV",
            data=csv,
            file_name=f"{exp_name}.csv",
            mime='text/csv',
        )
    else:
        st.write('Choose conditions and press "Count"')


    # c3, c4 = st.columns(2)
    # condition1 = c3.radio('Choose first condition',[f.name.split('.')[0] for f in ic.bed_files])
    # condition2 = c4.radio('Choose second condition', [f.name.split('.')[0] for f in ic.bed_files])
    # if condition1 == condition2:
    #     st.write('You chose the same condition')
    # else:
    #     df_list = list(ic.overlaps.values())
    #     st.write(df_list[0].head())
