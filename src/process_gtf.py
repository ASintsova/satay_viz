import pandas as pd

def change_chr_names(x):
    if x == 'Mito':
        return x
    else:
        return 'chr'+x


def process_gtf(gtf_file, out_file, feature="gene"):
    ann_df = pd.read_table(gtf_file, header=None, sep='\t', comment='#')
    ann_df = ann_df[ann_df[2] == feature]
    ann_df = ann_df[[0, 3, 4, 6, 8]]
    ann_df.columns = ['Chromosome', 'Start', 'End', 'Strand', 'gene_info']
    attributes = ["gene", "locus_tag"]
    for att in attributes:
        pattern = f'({att} .+?;|{att} .+?$)'
        ann_df[att] = (ann_df['gene_info'].str.extract(pattern, expand=False)
                                               .str.replace(f'{att} ', '')
                                               .str.strip(';""'))
    ann_df['gene'] = ann_df['gene'].fillna(ann_df['locus_tag'])
    ann_df['Chromosome'] = ann_df['Chromosome'].apply(change_chr_names)
    ann_df[['Chromosome', 'Start', 'End', 'Strand', 'gene', 'locus_tag']].to_csv(out_file, sep='\t', index=False)


if __name__ == "__main__":
    gtf_file = '/Users/ansintsova/git_repos/parfenova_satay/data/GCF_000146045.2_R64_genomic.edited.gtf'
    process_gtf(gtf_file, out_file="./annotations.csv", feature="gene")