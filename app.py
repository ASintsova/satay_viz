import streamlit as st

import plotly.graph_objects as go
#results_file = "/Users/ansintsova/git_repos/parfenova_satay/data/02_22_tn_counts/tn_counts_over_500_intervals_results.csv"
import pandas as pd
import plotly.express as px
import numpy as np

st.set_page_config(page_title="SATAY", layout="wide")

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


st.write("## Intervals close to gene of interest")

st.write('### Show intervals and there log2 fold changes near a gene of interest')
st.write("Each horizontal line shows corresponds to an interval and shows its log2 fold change" )

gene_list = st.multiselect('Choose gene(s) of interest', list(fres.gene_id.unique()), default=['YAL040C'])

if not gene_list:
    st.stop()

def extract_gene_info(df, gene_list, cutoff=0.05):
    gdf = df[df.gene_id.isin(gene_list)].copy()
    gdf['hits'] = gdf.padj < cutoff
    gdf['clrs'] = [px.colors.qualitative.Plotly[1] if c else px.colors.qualitative.Plotly[0] for c in gdf.hits.values]

    gene_info = {}
    for gene in gene_list:
        one_gene_df = gdf[gdf.gene_id == gene]
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


gdf, gene_info = extract_gene_info(fres, gene_list)
gene_def = draw_gene(gene_info)
line_info = gene_info['intervals']

x_text = []
for d in gene_def:
    if d['type'] == 'rect':
        x_text.append(d['x0']+ ((d['x1'] - d['x0'])/2))


fig = go.Figure()
# Set axes properties
fig.update_xaxes(range=[min(line_info['x0_intervals']) - 100, max(line_info['x1_intervals']) + 100], showgrid=False,
                 title='Position, bp')
fig.update_yaxes(range=[-3, 3], title='log2 fold change')

fig.update_layout(
    # filled Triangle
    height=600,
    shapes=gene_def,

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