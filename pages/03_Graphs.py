import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import numpy as np


def show_chrom_view(res_df, name, chrom, pval_th):
    logpval = f'{name}_logpval'
    hits = f'{name}_hits'
    foldChange = f'{name}_foldChange'
    df = res_df[res_df['Chromosome'] == chrom].copy()
    df[logpval] = -10 * np.log10(df.padj)
    df[hits] = df.padj < pval_th
    df = df.rename({'foldChange': foldChange}, axis=1)
    fig = px.scatter(df, x='Start', y=foldChange, size=logpval, hover_data={'location': False,
                                                                                'gene': True,
                                                                                'Distance': True,
                                                                                logpval: False},
                     hover_name='location', color=hits,
                     template='simple_white', labels={'Start': 'Position, bp', foldChange: f'{name} LFC',
                                                      'gene': 'closest gene', 'Distance': 'distance to gene',
                                                      hits: 'Hits'},
                     color_discrete_map={True: px.colors.qualitative.Plotly[1], False: px.colors.qualitative.Plotly[0]},
                     height=800, width=2000)

    fig.update_traces(marker=dict(
        opacity=0.7,
        line=dict(width=1,
                  color='DarkSlateGrey')),
        selector=dict(mode='markers'))
    fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, font=dict(size=18))
    return df, fig


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
            point=(one_gene_df.gene_start.values[0], -0.1), strand=list(one_gene_df.Strand.unique())[0])

    gene_info['intervals'] = dict(y_intervals=gdf.foldChange.values,
                                  x0_intervals=gdf.Start.values,
                                  x1_intervals=gdf.End.values, clrs=gdf.clrs.values)
    return gdf, gene_info


def draw_gene(gene_info):
    gdf_shape = []
    for gene, val in gene_info.items():
        if gene == 'intervals':
            continue
        if val['strand'] == '-':
            gdf_shape.append(dict(
                type="path",
                path=f" M {val['x0'] - val['rwidth'] / 10} 0 L {val['x0']} 0.2 L {val['x0']} -0.2 Z",
                fillcolor="LightGrey",
                line_color="LightGrey"))
        elif val['strand'] == '+':
            gdf_shape.append(dict(
                type="path",
                path=f" M {val['x1'] + val['rwidth'] / 10} 0 L {val['x1']} 0.2 L {val['x1']} -0.2 Z",
                fillcolor="LightGrey",
                line_color="LightGrey"))
        gdf_shape.append(dict(type="rect", x0=val['x0'], y0=-0.1, x1=val['x0'] + val['rwidth'], y1=0.1,
                              line=dict(color="LightGrey", width=2, ),
                              fillcolor="LightGrey", ))
    return gdf_shape


def draw_gene_and_intervals(line_info, gdf, gene_def, x_text, gene_list):
    layout = go.Layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )

    fig = go.Figure(layout=layout)
    # Set axes properties
    fig.update_xaxes(range=[min(line_info['x0_intervals']) - 100, max(line_info['x1_intervals']) + 100], showgrid=False,
                     title='Position, bp')
    fig.update_yaxes(range=[gdf.foldChange.min() - 0.5, gdf.foldChange.max() + 0.5], title='log2 fold change')
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
    return fig


def app():
    st.write("# Visualizing changes in #s of transposons over short intervals")
    st.markdown("#### Step 1: Upload the LFC/z-score file generated on `Analysis` page")
    results_file = st.file_uploader('Upload interval-results file')
    if not results_file:
        st.stop()
    st.markdown("#### Step 2: Enter the name/label for the comparison")
    exp_name = st.text_input('Enter comparison name', value=results_file.name)
    fres = pd.read_csv(results_file)
    cols_needed = ['Chromosome', 'Start', 'foldChange', 'location',
                   'gene', 'Distance', 'gene_start', 'gene_end', 'Strand']
    if not all([c in fres.columns for c in cols_needed]):
        st.write('Wrong file format')
        st.stop()

    with st.expander("Optional: Upload second comparision file and enter name/label"):
        results_file2 = st.file_uploader('Upload second interval-results file')
        if results_file2:
            exp_name2 = st.text_input('Enter second comparison name', value=results_file2.name)
            fres2 = pd.read_csv(results_file2)
            cols_needed = ['Chromosome', 'Start', 'foldChange', 'location',
                           'gene', 'Distance', 'gene_start', 'gene_end', 'Strand']
            if not all([c in fres2.columns for c in cols_needed]):
                st.write(f'Wrong file format for {results_file2.name}. Will be ignored')
                results_file2 = ''
            elif results_file.name == results_file2.name:
                st.write("Uploaded the same file twice. Will be ignored")
                results_file2 = ''

    st.markdown("""#### Chromosome-level View""")
    with st.expander("Chromosome-level View"):
        st.write(" Each interval is shown as a point, with genomic position on the x-axis, "
                 "and lfc between condition1 and conditon2 on the y-axis. "
                 "The size of the circle corresponds -10log(adjusted pvalue) "
                 " Hovering over each dot, you can see further information about the interval, "
                 "eg. the closest downstream gene and distance to that gene. "
                 "You can choose which chromosome to display, and the p-value cutoff that will define 'hits'.")

        c1, c2 = st.columns(2)
        chrom = c1.selectbox('Choose chromosome to display', list(fres['Chromosome'].unique()))
        pval_th = c2.number_input('Highlight intervals with p-values (padj) < ', value=0.05)
        df1, fig1 = show_chrom_view(fres, exp_name, chrom, pval_th)
        st.subheader(exp_name)
        st.plotly_chart(fig1, use_container_width=True)
        if results_file2:
            df2, fig2 = show_chrom_view(fres2, exp_name2, chrom, pval_th)
            st.subheader(exp_name2)
            st.plotly_chart(fig2, use_container_width=True)
            df3 = df1.merge(df2[['location', f'{exp_name2}_foldChange',
                                 f'{exp_name2}_logpval',
                                 f'{exp_name2}_hits']], how='inner', on='location')
            df3['hits2'] = df3[f'{exp_name}_hits'].astype(int) + df3[f'{exp_name2}_hits'].astype(int)
            hit_labels = {0: 'Not a hit', 1: 'Hit in one of the comparisons', 2: 'Hit in both comparisons'}
            df3['hits'] = df3['hits2'].map(hit_labels)
            st.subheader('Hits in only one of the 2 conditions:')
            st.write(df3[df3.hits == 'Hit in one of the comparisons'][['location', f'{exp_name}_foldChange',
                                                                       f"{exp_name2}_foldChange",
                                                                       'gene', 'locus_tag', 'Distance']])
            fig3 = px.scatter(df3, x=f"{exp_name}_foldChange", y=f'{exp_name2}_foldChange',
                              color='hits', height=800, width=800, hover_data=['location', 'gene', 'locus_tag',  'Distance'],
                              color_discrete_map={'Not a hit': 'grey', 'Hit in one of the comparisons': px.colors.qualitative.Plotly[1],
                                                  'Hit in both comparisons': px.colors.qualitative.Plotly[0]})
            fig3.update_traces(marker=dict(
                opacity=0.7,
                size=12,
                line=dict(width=1,
                          color='DarkSlateGrey')),
                selector=dict(mode='markers'))
            fig3.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, font=dict(size=18))
            fig3.update_xaxes(range=[min(df3[f'{exp_name}_foldChange'].min(), df3[f'{exp_name2}_foldChange'].min()) - 0.5,
                                     max(df3[f'{exp_name}_foldChange'].max(), df3[f'{exp_name2}_foldChange'].max()) + 0.5])
            fig3.update_yaxes(range=[min(df3[f'{exp_name}_foldChange'].min(), df3[f'{exp_name2}_foldChange'].min()) - 0.5,
                                     max(df3[f'{exp_name}_foldChange'].max(), df3[f'{exp_name2}_foldChange'].max()) + 0.5])
            fig3.add_hline(y=0, line_dash="dash", line_color="grey")
            fig3.add_vline(x=0, line_dash="dash", line_color="grey")
            st.write(f"This graph shows LFC in {exp_name} vs LFC in {exp_name2}.")
            st.plotly_chart(fig3, use_container_width=True)


    with st.expander("Gene-level View"):
        st.markdown("This graph shows intervals (shown as horizontal lines) and their lfcs near a gene of interest.")
        gene_list = st.multiselect('Choose gene(s) of interest', list(fres.gene.unique()))
        if not gene_list:
            st.stop()

        gdf, gene_info = extract_gene_info(fres, gene_list, cutoff=pval_th)
        gene_def = draw_gene(gene_info)
        line_info = gene_info['intervals']

        x_text = []
        for d in gene_def:
            if d['type'] == 'rect':
                x_text.append(d['x0'] + ((d['x1'] - d['x0'])/2))

        fig4 = draw_gene_and_intervals(line_info, gdf, gene_def, x_text, gene_list)
        st.subheader(exp_name)
        st.plotly_chart(fig4, use_container_width=True)
        c3, c4 = st.columns(2)
        if results_file2:
            if all([gene in list(fres2.gene.unique()) for gene in gene_list]):
                gdf2, gene_info2 = extract_gene_info(fres2, gene_list, cutoff=pval_th)
                gene_def2 = draw_gene(gene_info2)
                line_info2 = gene_info2['intervals']
                x_text2 = []
                for d in gene_def2:
                    if d['type'] == 'rect':
                        x_text2.append(d['x0'] + ((d['x1'] - d['x0']) / 2))

                fig5 = draw_gene_and_intervals(line_info2, gdf2, gene_def2, x_text2, gene_list)
                st.subheader(exp_name2)
                st.plotly_chart(fig5, use_container_width=True)
            else:
                c4.write(f"{gene_list} not found in {results_file2.name}")

        c3.write(exp_name)

        c3.write(gdf.sort_values('Start')[['location', 'foldChange',
                                                    'padj', 'gene', 'Distance', 'gene_start',
                                                    'gene_end', 'Strand']])
        if results_file2 and all([gene in list(fres2.gene.unique()) for gene in gene_list]):
            c4.write(exp_name2)
            c4.write(gdf2.sort_values('Start')[['location', 'foldChange',
                                            'padj', 'gene', 'Distance', 'gene_start',
                                            'gene_end', 'Strand']])

app()
