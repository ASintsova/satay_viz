from src.count_inserts import IntervalCounter
from src.compare_interval_counts import analyze_one_offset, get_tts, concatenate_results
import streamlit as st
import pandas as pd
from datetime import datetime as dt
import pyranges as pr


@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')


def app():
    st.header("Compare transposon insertion frequency between two conditions")
    today = dt.now().strftime("%d-%m-%Y")
    st.markdown("#### Step 1: Upload the transposon count file (generated on `Count` page)")
    input_file = st.file_uploader('Upload transposon counts', accept_multiple_files=False)
    st.markdown("#### Step 2: Choose conditions to compare")
    if not input_file:
        st.stop()
    count_df = pd.read_csv(input_file, index_col=0)
    fixed_col = ['Chromosome', 'Start', 'End', 'offset']
    samples = [c for c in count_df.columns if c not in fixed_col]
    c1, c2 = st.columns(2)
    condition1 = c1.radio('Choose first condition', samples)
    condition2 = c2.radio('Choose second condition', [s for s in samples if s != condition1])
    st.markdown("#### Step 3: Calculate LFC and z-scores for all intervals.")
    if st.button('Compare'):
        st.write(f'Calculating LFC for {condition1} vs {condition2} ...')
        df_list = []
        for offset, df in count_df.groupby('offset'):
            _, results = analyze_one_offset(df, offset, condition1, condition2, get_tts)
            df_list.append(results)
        fdf = pd.concat(df_list)
        fdf = fdf.reset_index()
        fdf['Chromosome'] = fdf.location.str.split(":", expand=True)[0]
        positions = fdf.location.str.split(":", expand=True)[1].str.split('-', expand=True)
        fdf['Start'] = positions[0].astype(int)
        fdf['End'] = positions[1].astype(int)
        intervals = pr.PyRanges(fdf)
        st.write('Annotating results...')
        genes = pr.PyRanges(pd.read_table('annotations.tsv'))
        nearest = intervals.nearest(genes).as_df().rename({'Start_b': 'gene_start', 'End_b': 'gene_end'}, axis=1)
        st.write("Done")
        csv = convert_df(nearest)
        st.markdown("#### Step 4: Download and save the results.")
        st.download_button(
            label="Download counts as CSV",
            data=csv,
            file_name=f"{input_file.name}_{condition1}_vs_{condition2}.csv",
            mime='text/csv',
        )
    else:
        st.write('Choose conditions and press "Compare"')

app()
