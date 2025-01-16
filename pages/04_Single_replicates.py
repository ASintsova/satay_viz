import streamlit as st
import pandas as pd
from upsetplot import plot
import matplotlib.pyplot as plt
import numpy as np


def generate_upset_plot(df, cat_order):
    unique_values = pd.unique(df.values.ravel())
    unique_values = [v for v in unique_values if v is not np.nan]
    result = pd.DataFrame(
        {val: df.apply(lambda col: col == val).any(axis=0) for val in unique_values},
        index=df.columns
    ).T
    result = result[cat_order]
    result_vc = result.value_counts()
    plot(result_vc, sort_categories_by='-input')
    fig = plt.gcf()  # Get the current figure
    fig.set_dpi(300)
    fig.set_size_inches(6, 4) 
    return fig, result


def app():
    st.header("Visualizing single replicate analysis")
    st.markdown("## UpSet Plots")
    st.markdown("#### Upload hits for each replicate in the following format: ")
    upset_df = pd.DataFrame([
                               "a,b,b,d".split(","), "b,a,c,e".split(","), "f,g,f,h".split(","),
                               "i,i,i,i".split(","), "g,k,,l".split(","), "m,n,,".split(",")], 
                               columns=['Rep1', 'Rep2', 'Rep3', 'Pooled'])
    st.markdown("The file can have any number of columns, and column labels do NOT need to match the example")
    st.markdown("Here is an example dataframe:")
    st.dataframe(upset_df)
    hits_file = st.file_uploader('File upload')
    if not hits_file:
        st.warning("No file uploaded. Showing an example UpSet plot.")
    else:
        upset_df = pd.read_csv(hits_file)
    selection = st.pills("Select category order:", upset_df.columns, selection_mode="multi", default=upset_df.columns.sort_values())
    st.markdown(f"Your selected category orders: {selection}.")
    if len(selection) < 2:
        st.warning("Cannot draw upset plot for 1 category, select more categories")
        st.stop()
    fig, result = generate_upset_plot(upset_df, selection)
    st.pyplot(fig)
    st.markdown("#### Show specific intersections")
    reps = st.multiselect(label='Choose intersection to show', options=['All'] + list(result.columns), default='All')
    if 'All' in reps:
        st.write(result[result.index != ''])
    else:
        # Filter rows
        filtered = result[
            result[reps].all(axis=1) &  # All specified columns are True
            result.drop(columns=reps).eq(False).all(axis=1)  # All other columns are False
        ]
        st.write(filtered)
app()