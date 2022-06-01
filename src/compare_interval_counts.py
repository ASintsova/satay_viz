import pandas as pd
import numpy as np
import scipy
from scipy.stats import zscore
from statsmodels.stats import multitest
#from skbio.stats.composition import clr
from pathlib import Path

def get_tts(df):
    return np.log2(df/df.sum()*1000000 + 0.5)

def calculate_zscore(cnts, col1, col2):
    results = cnts.copy()
    results['foldChange'] = cnts[col1] - cnts[col2]
    results['mean_cnts'] = results[[col1, col2]].mean(axis=1)
    results['zscore'] = zscore(results.foldChange)
    results['pval'] = results.zscore.apply(lambda x: scipy.stats.norm.sf(abs(x)) * 2)
    results['padj'] = multitest.multipletests(results['pval'], alpha=0.05, method='fdr_bh')[1]
    return results


def analyze_one_offset(count_df, offset, condition1, condition2, norm_func):
    df = count_df.copy()
    df['location'] = df['Chromosome'] + ":" + df['Start'].astype(str) + "-" + df['End'].astype(str)
    df = df[['location', condition1, condition2]].set_index('location')
    df = df[df.sum(axis=1) > 10]
    norm_df = norm_func(df)
    results = calculate_zscore(norm_df, condition1, condition2).assign(offset=offset)
    return df, results


def concatenate_results(df_list, condition1, condition2, norm_func=get_tts):
    results_list = []
    for df in df_list:
        _, results = analyze_one_offset(df, condition1, condition2, norm_func)
        results_list.append(results)
    fdf = pd.concat(results_list)
    return fdf
    #fdf.to_csv(Path(outdir)/f"{condition1}_vs_{condition2}_results.csv")

