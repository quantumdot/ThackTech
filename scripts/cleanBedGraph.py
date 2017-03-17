#!/usr/bin/env python
import os
import sys
import argparse
from ThackTech import bdgtools, chromtools



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bedgraph', help='BedGraph file to clean.')
    parser.add_argument('-o', '--output', help='Specifies the output file. Default is stdout')
    parser.add_argument('-f', '--format', choices=['bdg', 'bw'], default='bdg', help='Output format.')
    parser.add_argument('-g', '--genome', action='store', required=True, help='Genome chromosome sizes. Can specify UCSC genome builds (i.e. hg19, mm9) or the location of a chromosome sizes file (standard UCSC genome sizes format).')
    parser.add_argument('-r', '--repairoverlaps', action='store_true', help='Repair overlapping intervals in the bedgraph, substituting the intersecting regions with a new interval having a score computed by --method.')
    parser.add_argument('-m', '--method', choices=bdgtools.score_funcs.keys(), default='mean', help='Method to use for computation of score when repairing/merging overlapping intervals.')
    parser.add_argument('-e', '--missingregions', action='store', default='ignore', help='How to treat regions with no data. `ignore` to ignore missing regions, otherwise a numeric (float or int).')
    parser.add_argument('-q', '--quiet', action='store_true', help='Be as quiet as possible.')
    args = parser.parse_args()
    
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
        
    if args.quiet:
        sys.stderr = os.devnull
        
    chrom_sizes = chromtools.ChromSizes(args.genome)
    #chrom_sizes = bdgtools.get_chrom_sizes(args.genome)
        

    with open(args.bed, 'r') as infile:
        intervals = bdgtools.parse_bedgraph(infile)
    
    
    if args.repairoverlaps:
        intervals = bdgtools.reduce_overlaps(intervals, chrom_sizes, bdgtools.score_funcs[args.method])
    
    
    if not args.missingregions == 'ignore':
        intervals = bdgtools.fill_complement(intervals, chrom_sizes, float(args.missingregions))

    
    if args.format == 'bw':
        bdgtools.write_bigwig(intervals, chrom_sizes, sys.stdout)
    else:
        bdgtools.write_bedgraph(intervals, sys.stdout)
    return
#end main()


if __name__ == "__main__":
    main()
