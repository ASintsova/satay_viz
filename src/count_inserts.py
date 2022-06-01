import pyranges as pr
import pandas as pd
import typing
from typing import Union, List, Dict
from pathlib import Path
import sys


class IntervalCounter:

    def __init__(self, bed_files: List, chr_sizes_file: str, interval_size: int, interval_offset: int,
                 name: str, outdir: str = ".", sep: str = " ", skiprows: int = 1):
        self.bed_files = bed_files
        self.chr_sizes_file = chr_sizes_file
        self.chr_sizes = {}
        self.name = name
        self.sep = sep
        self.skiprows = skiprows
        self.outdir = Path(outdir)
        self.interval_size = interval_size
        self.interval_offset = interval_offset
        self.offsets = list(range(0, self.interval_size, self.interval_offset))
        self.custom_ranges = []
        self.tn_sites = {}
        self.overlaps = {}
        self.out_files = [self.outdir/f"{self.name}_offset{s}.csv" for s in self.offsets]

    def read_chr_sizes(self):
        chr_sizes_df = pd.read_table(self.chr_sizes_file)
        if 'Chromosome' not in chr_sizes_df.columns or 'Length' not in chr_sizes_df.columns:
            print('Cannot process Chromosome sizes file. Make sure columns "Chromosome" and "Length" are present')
            sys.exit(1)
        chr_names = chr_sizes_df.Chromosome.values
        chr_sizes = chr_sizes_df.Length.values
        self.chr_sizes = {name: int(size) for name, size in zip(chr_names, chr_sizes)}

    # def generate_chr_range(self, chr_name: str, chr_len: int, start: int = 0, size: int = 500) -> dict:
    #     """
    #     Given a chromosome name and length, generate ranges of a given size, starting at a given position
    #     :param chr_name: chromosome name
    #     :param chr_len: chromosome length
    #     :param start: starting position for the intervals
    #     :param size: size of the intervals
    #     :return: dictionary of intervals
    #     """
    #     starts = list(range(start, chr_len, size))
    #     ends = list(range(start + size, chr_len, size)) + [chr_len]
    #     chrs = [chr_name] * len(starts)
    #     return {'Chromosome': chrs, "Start": starts, "End": ends}

    def create_custom_ranges(self):
        """
        Given all a dict of chromosome names and length, create ranges of given size for each, concat into a df
        :param chr_map: dict of {name: length}
        :param start: starting position for the intervals
        :param size: size of the intervals
        :return: ranges dataframe
        """
        for offset in self.offsets:
            genomic_ranges = []
            for chr_name, chr_len in self.chr_sizes.items():
                interval_starts = list(range(offset, chr_len, self.interval_size))
                interval_ends = list(range(offset + self.interval_size, chr_len, self.interval_size)) + [chr_len]
                chrs = [chr_name] * len(interval_starts)
                chr_ranges = pr.from_dict({'Chromosome': chrs, "Start": interval_starts, "End": interval_ends})
                genomic_ranges.append(chr_ranges)
            ranges_gr = pr.concat(genomic_ranges)
            self.custom_ranges.append(ranges_gr)

    def get_tn_sites_for_all_samples(self):
        """
        Take list of sample files and convert to dictionary of genomic ranges
        :param sample_files:
        :return: dictionary of tn insertion sites for each sample
        """
        col_names = "Chromosome Start End Strand Score".split()
        for bed_file in self.bed_files:
            name = bed_file.name.replace("-", "_").split('.')[0]
            print(name)
            df = pd.read_table(bed_file, sep=self.sep, skiprows=self.skiprows)
            df.columns = col_names[0:df.shape[1]]
            self.tn_sites[name] = pr.PyRanges(df)

    def find_overlaps(self):
        """
        Take dictionary of tn insertion sites and genomic ranges and calculate overlaps
        :param tn_sites:
        :param chr_gr:
        :return:
        """
        for custom_range, out_file in zip(self.custom_ranges, self.out_files):
            overlap = pr.count_overlaps(self.tn_sites, custom_range)
            #overlap.to_csv(out_file)
            self.overlaps[out_file.name] = overlap

    def count(self):
        # Get chromosome sizes
        self.read_chr_sizes()
        # Generate custom ranges for each offset
        self.create_custom_ranges()
        # Get tn sites for all samples
        self.get_tn_sites_for_all_samples()
        # Find overlaps between custom range and tn sites
        self.find_overlaps()

"""

Part 2: Count reads over different intervals, generate 9 count files

1. Generate gtf/saf file of size 500 across genome

GeneID	Chr	Start	End	Strand
497097	chr1	3204563	3207049	-
497097	chr1	3411783	3411982	-
497097	chr1	3660633	3661579	-


2. featureCounts the reads
3. ouptut count file
4. Do this 9 times

"""
# if __name__ == "__main__":
#     data_dir = Path("/data/")
#     out_dir = data_dir/"02_22_tn_counts"
#     samples = [data_dir/"15776-1_noaF/15776-1_noaF.bam.bed", data_dir/"15776-1_4nMaF/15776-1_4nMaF.bam.bed"]
#     starts = [1] + list(range(50, 500, 50))
#     size = 500
#     for start in starts:
#         overlaps = count_tns_over_custom_range(chr_map, start, size, samples, out_dir/f"15776-1-4aF-0aF-{size}-{start}.csv")

# chr_names = """Chromosome	Genbank ID	RefSeq ID	Length (bp)
# # Chromosome I	BK006935.2	NC_001133.9	230218
# # Chromosome II	BK006936.2	NC_001134.8	813184
# # Chromosome III	BK006937.2	NC_001135.5	316620
# # Chromosome IV	BK006938.2	NC_001136.10	1531933
# # Chromosome V	BK006939.2	NC_001137.3	576874
# # Chromosome VI	BK006940.2	NC_001138.5	270161
# # Chromosome VII	BK006941.2	NC_001139.9	1090940
# # Chromosome VIII	BK006934.2	NC_001140.6	562643
# # Chromosome IX	BK006942.2	NC_001141.2	439888
# # Chromosome X	BK006943.2	NC_001142.9	745751
# # Chromosome XI	BK006944.2	NC_001143.9	666816
# # Chromosome XII	BK006945.2	NC_001144.5	1078177
# # Chromosome XIII	BK006946.2	NC_001145.3	924431
# # Chromosome XIV	BK006947.3	NC_001146.8	784333
# # Chromosome XV	BK006948.2	NC_001147.6	1091291
# # Chromosome XVI	BK006949.2	NC_001148.4	948066
# # Chromosome Mito	AJ011856.1	NC_001224.1	85779""".split('\n')
# #
# # chr_names = [ch.split("\t") for ch in chr_names if 'Mito' not in ch]
# # chr_map = {f"chr{c[0].split()[1]}": int(c[3]) for c in chr_names[1:]}