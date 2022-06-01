import click
import subprocess
import shlex
from pathlib import Path
import click
import sys

from src.count_inserts import IntervalCounter
from src.compare_interval_counts import concatenate_results

class DefaultHelp(click.Command):
    def __init__(self, *args, **kwargs):
        context_settings = kwargs.setdefault('context_settings', {})
        if 'help_option_names' not in context_settings:
            context_settings['help_option_names'] = ['-h', '--help']
        self.help_flag = context_settings['help_option_names'][0]
        super(DefaultHelp, self).__init__(*args, **kwargs)

    def parse_args(self, ctx, args):
        if not args:
            args = [self.help_flag]
        return super(DefaultHelp, self).parse_args(ctx, args)

@click.group(options_metavar='', subcommand_metavar='<command> [options]', invoke_without_command=False)
def main():
    """
    \b
    Program: satay - a tool for analysis of satay data over custom intervals
    Version: 1.0.0
    Reference:

    Type satay <command> to print the help for a specific command

    """

####################
#   Count Inserts  #
####################

@main.command(cls=DefaultHelp, short_help="\tcount transposons over custom intervals", options_metavar='<options>')
@click.option('--input_files', '-i', default='', help='input bed file specifying positin of each transposon in a sample.', metavar='FILE[,FILE]')
@click.option('--bed_dir', '-d', default='',  help='process all bed files in the directory', metavar='DIR')
@click.option('--name', '-n', default='Test', help='Output file prefix', metavar='STR')
@click.option('--chr_sizes', '-c', required=False, default='chr_sizes.tsv', help='TSV file specifying chromosome names and sizes', metavar='FILE')
@click.option('--interval_size', '-s', default=500, help='interval size', metavar='INT')
@click.option('--offset', '-f', default=50, help='offset between intervals', metavar='INT')
@click.option('--out_dir', '-o', default='.', help='output directory [.]', metavar="DIR")
@click.option('--sep', '-t', default=' ', help='separator to use for bed files', metavar="STR")
@click.option('--skiprows', '-l', default=1, help='number of header rows in the bed files', metavar="INT")
def count(input_files, bed_dir, chr_sizes,  name, interval_size, offset,  out_dir, sep, skiprows):
    if not input_files and not bed_dir:
        sys.exit('Please specify bed files to process or a directory with bed files')
    elif input_files:
        input_files = [Path(f) for f in input_files.split(',')]
    else:
        input_files = [f for f in Path(bed_dir).glob("**/*.bed")]

    intDS = IntervalCounter(input_files, chr_sizes, interval_size, offset, name, out_dir,
                            sep, skiprows)
    intDS.count()


@main.command(cls=DefaultHelp, short_help="\tcompare transposons counts between 2 conditions", options_metavar='<options>')
@click.option('--input_files', '-i', default='', help='input offset files produced wtih count.', metavar='FILE[,FILE]')
@click.option('--file_dir', '-d', default='',  help='process all offset files in the directory', metavar='DIR')
@click.option('--condition1', '-s1', help='first condition to compare', metavar='STR')
@click.option('--condition2', '-s2',  help='second condition to compare', metavar='INT')
@click.option('--out_dir', '-o', default='.', help='output directory [.]', metavar="DIR")
def compare(input_files, file_dir, condition1, condition2, out_dir):
    if not input_files and not file_dir:
        sys.exit('Please specify bed files to process or a directory with bed files')
    elif input_files:
        input_files = [Path(f) for f in input_files.split(',')]
    else:
        input_files = [f for f in Path(file_dir).glob("**/*offset*.csv")]
    concatenate_results(input_files, condition1, condition2, out_dir)

