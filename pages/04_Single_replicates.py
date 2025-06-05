import streamlit as st
import pandas as pd
from upsetplot import plot
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px

def get_annotations(annotation_file="annotations.tsv"):
    return pd.read_table(annotation_file)

def show_ranking_plot(res_df, metric_type):
    """
    Show genes ranked from most negative to most positive based on selected metric
    """
    df = res_df.dropna(subset=[metric_type]).copy()
    df = df.sort_values(metric_type)
    df['rank'] = range(len(df))
    
    fig = px.scatter(df, x='rank', y=metric_type, 
                     hover_data={'locus_tag': False, 'gene': True, 'rank': False},
                     hover_name='locus_tag', color='hits',
                     template='simple_white', 
                     labels={'rank': 'Gene Rank (most negative to positive)', 
                             'foldChange': 'Fold Change',
                             'zscore': 'Z-score',
                             'hits': 'Hits'},
                     color_discrete_map={True: px.colors.qualitative.Plotly[1], False: px.colors.qualitative.Plotly[0]},
                     height=500, width=1200)
    
    fig.update_traces(marker=dict(
        size=4,
        opacity=0.7,
        line=dict(width=1, color='DarkSlateGrey')),
        selector=dict(mode='markers'))
    fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, font=dict(size=18))
    return fig

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


def show_chrom_view(res_df,  chrom,):
    """
    # TODO repetitive with Intervals.py
    """
    df = res_df[res_df['Chromosome'] == chrom].dropna(subset=['logpval'])
    fig = px.scatter(df, x='Start', y='foldChange', size='logpval', hover_data={'locus_tag': False,
                                                                                'gene': True,
                                                                            
                                                                                'logpval': False},
                     hover_name='locus_tag', color='hits',
                     template='simple_white', labels={'Start': 'Position, bp', 'log2FoldChange': f'LFC',
                                                      'gene': 'closest gene', 'Distance': 'distance to gene',
                                                      'hits': 'Hits'},
                     color_discrete_map={True: px.colors.qualitative.Plotly[1], False: px.colors.qualitative.Plotly[0]},
                     height=500, width=2000)

    fig.update_traces(marker=dict(
        size=4,
        opacity=0.7,
        line=dict(width=1,
                  color='DarkSlateGrey')),
        selector=dict(mode='markers'))
    fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, font=dict(size=18))
    return df, fig



def app():

    gene_ids = get_annotations()
    st.header("Visualizing single replicate analysis")

    st.markdown('## Gene LFC')

    uploaded_files = st.file_uploader("Upload multiple files", accept_multiple_files=True)
    file_dfs = {}
    
    # User controls for defining hits
    st.sidebar.header("Hit Definition")
    threshold_type = st.sidebar.radio("Define hits based on:", ["foldChange", "zscore"])
    threshold_value = st.sidebar.number_input(f"{threshold_type} threshold (absolute value)", min_value=0.0, value=2.0, step=0.1)
    
    if uploaded_files:
        for uploaded_file in uploaded_files:
            df = pd.read_csv(uploaded_file, index_col=0)
            # Add analysis column if missing
            if 'analysis' not in df.columns:
                if 'tn' in uploaded_file.name.lower():
                    df['analysis'] = 'tn'
                elif 'read' in uploaded_file.name.lower():
                    df['analysis'] = 'read'
                else:
                    df['analysis'] = 'unknown'
            df['logpval'] = 6  # Set logpval to 6 for all
            df['hits'] = abs(df[threshold_type]) > threshold_value
            merge_on = [c for c in df.columns if c in gene_ids.columns]
            file_dfs[uploaded_file.name] = df.merge(gene_ids, on=merge_on, how='left')
        st.write(file_dfs[uploaded_file.name].head())
    else:
        # Load example file if no files uploaded
        st.info("No files uploaded. Loading example file: Rep1_4_vs_0nMaF_tn_subsample.csv.gz")
        example_df = pd.read_csv("examples/Rep1_4_vs_0nMaF_tn_subsample.csv.gz", index_col=0)
        # Add analysis column if missing
        if 'analysis' not in example_df.columns:
            example_df['analysis'] = 'tn'  # Example file has 'tn' in name
        example_df['logpval'] = 6  # Set logpval to 12 for all
        example_df['hits'] = abs(example_df[threshold_type]) > threshold_value
        merge_on = [c for c in example_df.columns if c in gene_ids.columns]
        file_dfs["Rep1_4_vs_0nMaF_tn.csv"] = example_df.merge(gene_ids, on=merge_on, how='left')
        st.write(file_dfs["Rep1_4_vs_0nMaF_tn.csv"].head())
            
    comps = file_dfs.keys()
    
    # Get all unique analysis types
    all_analysis_types = set()
    for df in file_dfs.values():
        all_analysis_types.update(df['analysis'].unique())
    
    # Analysis type filter
    if len(all_analysis_types) > 1:
        analysis_filter = st.selectbox('Filter by analysis type:', ['All'] + sorted(list(all_analysis_types)))
    else:
        analysis_filter = 'All'

    to_show = st.radio('Which comp to show', comps)
    chr_to_show = st.selectbox('Chromosome', gene_ids.Chromosome.unique())
    
    # Apply analysis filter to selected dataframe
    selected_df = file_dfs[to_show].copy()
    if analysis_filter != 'All':
        selected_df = selected_df[selected_df['analysis'] == analysis_filter]

    df, fig = show_chrom_view(selected_df, chr_to_show)

    st.plotly_chart(fig, use_container_width=True)
    
    # Add ranking plot
    st.markdown(f"## Gene Ranking by {threshold_type}")
    ranking_fig = show_ranking_plot(selected_df, threshold_type)
    st.plotly_chart(ranking_fig, use_container_width=True)
    
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