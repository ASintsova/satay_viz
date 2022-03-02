import streamlit as st

import plotly.graph_objects as go

import pandas as pd
import plotly.express as px
import numpy as np

st.set_page_config(page_title="SATAY", layout="wide")

st.write("# Visualizing changes in #s of transposons over 500-bp Intervals")
results_file = st.file_uploader('Upload interval-results file')
#results_file = "/Users/ansintsova/git_repos/parfenova_satay/data/02_22_tn_counts/tn_counts_over_500_intervals_results.csv"
if not results_file:
    st.stop()

fres = pd.read_csv(results_file)
cols_needed = ['chromosome', 'start', 'lfc', 'location', 'gene', 'distance']
if not all([c in fres.columns for c in cols_needed]):
    st.write('Wrong file format')
    st.stop()

st.markdown("""### Analysis Summary:
- Divide genome into 500 bp intervals (ex. chrI:1-500, chrI:501-100, etc.)
- Calculate # of transposons for each 500 nt interval for **15776-1-4nMaF** and **15776-1-noaF** samples.
- Normalize counts and calculate log2 Fold Change (lfc) between **15776-1-4nMaF** and **15776-1-noaF**
- Calclate z-score, pval and adjusted pval for each interval
- Repeat for intervals offset by 50 nt (ex. chrI:50-500, chrI:551:1050)

""")

st.markdown("""### Chromosome-level View:
Each interval is shown as a point, with genomic position on the x-axis, and lfc between 15776-1-4nMaF and 15776-1-noaF on the y-axis. The size of the circle corresponds -10log(adjusted pvalue).
Hovering over each dot, you can see further information about the interval, eg. the closest downstream gene and distance to that gene.
You can choose which chromosome to display, and the p-value cutoff that will define 'hits'.

""")

c1, c2 = st.columns(2)
chrom = c1.selectbox('Choose chromosome to display', list(fres['chromosome'].unique()))
pval_th = c2.number_input('Highlight intervals with p-values (padj) < ', value=0.05)
df = fres[fres['chromosome'] == chrom].copy()
df['logpval'] = -10*np.log10(df.padj)
df['hits'] = df.padj < pval_th
fig = px.scatter(df, x='start', y='lfc', size='logpval', hover_data={'location': False,
                                                                      'gene': True,
                                                                      'distance': True,
                                                                      'logpval':False},
                 hover_name='location', color='hits',
                 template='simple_white', labels={'start': 'Position, bp', 'lfc': 'log2 fold change',
                                           'gene': 'closest gene', 'distance': 'distance to gene'},
                 color_discrete_map={True: px.colors.qualitative.Plotly[1], False: px.colors.qualitative.Plotly[0]},
                 height=800, width=2000)

fig.update_traces(marker=dict(
                              opacity = 0.7,
                              line=dict(width=1,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'))
fig.update_layout({'paper_bgcolor':'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, font=dict(size=18))
st.plotly_chart(fig, use_container_width=True)


st.markdown("""### Gene-level View
This graph shows intervals (shown as horizontal lines) and their lfcs near a gene of interest. 
""")


gene_list = st.multiselect('Choose gene(s) of interest', list(fres.gene.unique()), default=['BUD14'])

if not gene_list:
    st.stop()


def extract_gene_info(df, gene_list, gene_col='gene', cutoff=0.05):
    gdf = df[df[gene_col].isin(gene_list)].copy()
    gdf['hits'] = gdf.padj < cutoff
    gdf['clrs'] = [px.colors.qualitative.Plotly[1] if c else px.colors.qualitative.Plotly[0] for c in gdf.hits.values]
    gene_info = {}
    for gene in gene_list:
        one_gene_df = gdf[gdf[gene_col] == gene]
        gene_info[gene] = dict(
            x0=one_gene_df.gene_start.values[0],
            x1=one_gene_df.gene_end.values[0],
            rwidth=abs(one_gene_df.gene_start.values[0] - one_gene_df.gene_end.values[0]),
            rheight=0.2,
            point=(one_gene_df.gene_start.values[0], -0.1), strand=list(one_gene_df.strand.unique())[0])

    gene_info['intervals'] = dict(y_intervals=gdf.lfc.values,
                                  x0_intervals=gdf.start.values,
                                  x1_intervals=gdf.end.values, clrs=gdf.clrs.values)
    return gdf, gene_info


def draw_gene(gene_info):
    gdf_shape = []
    for gene, val in gene_info.items():
        if gene == 'intervals':
            continue
        if val['strand'] == '-':
            gdf_shape.append(dict(
                type="path",
                path=f" M {val['x0']-val['rwidth']/10} 0 L {val['x0']} 0.2 L {val['x0']} -0.2 Z",
                fillcolor="LightGrey",
                line_color="LightGrey"))
        elif val['strand'] == '+':
            gdf_shape.append(dict(
                type="path",
                path=f" M {val['x1']+val['rwidth']/10} 0 L {val['x1']} 0.2 L {val['x1']} -0.2 Z",
                fillcolor="LightGrey",
                line_color="LightGrey"))
        gdf_shape.append(dict(type="rect", x0=val['x0'], y0=-0.1, x1=val['x0']+val['rwidth'], y1=0.1,
                              line=dict( color="LightGrey", width=2,),
                fillcolor="LightGrey",))
    return gdf_shape


gdf, gene_info = extract_gene_info(fres, gene_list, cutoff=pval_th)
gene_def = draw_gene(gene_info)
line_info = gene_info['intervals']

x_text = []
for d in gene_def:
    if d['type'] == 'rect':
        x_text.append(d['x0'] + ((d['x1'] - d['x0'])/2))


layout = go.Layout(
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)'
)
fig = go.Figure(layout=layout)

# Set axes properties
fig.update_xaxes(range=[min(line_info['x0_intervals']) - 100, max(line_info['x1_intervals']) + 100], showgrid=False,
                 title='Position, bp')
fig.update_yaxes(range=[gdf.lfc.min()-0.5, gdf.lfc.max()+0.5], title='log2 fold change')

fig.update_layout(
    # filled Triangle
    height=600,
    shapes=gene_def,
    template='plotly_white'

)

fig.add_trace(go.Scatter(
    x=x_text,
    y=[-0.25] * len(x_text),
    text=gene_list,
    mode="text",

))

for y, x0, x1, c in zip(line_info['y_intervals'], line_info['x0_intervals'], line_info['x1_intervals'],
                        line_info['clrs']):
    fig.add_shape(type="line", x0=x0,
                  y0=y,
                  x1=x1,
                  y1=y, line=dict(color=c), opacity=0.5)

st.plotly_chart(fig, use_container_width=True)

st.write(gdf.sort_values('start')[['location', 'lfc', 'padj', 'gene', 'distance', 'gene_start', 'gene_end', 'strand']])