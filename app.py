import streamlit as st

#results_file = "/Users/ansintsova/git_repos/parfenova_satay/data/02_22_tn_counts/tn_counts_over_500_intervals_results.csv"
import pandas as pd
import plotly.express as px
import numpy as np

st.set_page_config(page_title="Test",layout="wide")

st.write("# Visualizing Transposon Density over 500-bp Intervals")
results_file = st.file_uploader('Upload interval-results file')
if not results_file:
    st.stop()

fres = pd.read_csv(results_file)
cols_needed = ['chr2', 'start', 'lfc', 'location', 'gene_id', 'distance']
if not all([c in fres.columns for c in cols_needed]):
    st.write('Wrong file format')
    st.stop()

st.write('### Protocol:')
st.write('1. Calculate # of transposons for each 500 nt interval')
st.write('2. Normalize transposon counts')
st.write('3. Calculate log2 Fold Change between **15776-1-4nMaF** and **15776-1-noaF**')
st.write('4. Caclulate z-score for each lfc')


st.write('**For each interval, the plot shows the intervals position on the x-axis, and log2FC between 15776-1-4nMaF and 15776-1-noaF on the y-axis. The size of the circle corresponds -10log(adjusted pvalue)**')
st.write('Hovering over each dot, you can see further information about the interval, eg. the closest gene and distance to that gene')


c1, c2 = st.columns(2)
chrom = c1.selectbox('Choose chromosome to display', list(fres.chr2.unique()))
pval_th = c2.number_input('Highlight intervals with p-values (padj) < ', value=0.05)
df = fres[fres.chr2 == chrom].copy()
df['logpval'] = -10*np.log10(df.padj)
df['hits'] = df.padj < pval_th
fig = px.scatter(df, x='start', y='lfc', size='logpval', hover_data={'location': False,
                                                                       'gene_id': True,
                                                                       'distance': True,
                                                                      'logpval':False},
                 hover_name='location', color='hits',
                 template='simple_white', labels={'start':'Position, bp', 'lfc': 'log2 fold change',
                                           'gene_id': 'closest gene', 'distance': 'distance to gene'},
                 color_discrete_map={True: px.colors.qualitative.Plotly[1], False: px.colors.qualitative.Plotly[0]},
                 height=800, width=2000)

fig.update_traces(marker=dict(
                              opacity = 0.7,
                              line=dict(width=1,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'))
fig.update_layout(font=dict(size=18))
st.plotly_chart(fig, use_container_width=True)