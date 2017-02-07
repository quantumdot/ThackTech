#!/usr/bin/env python
import os
import sys
import argparse
import cStringIO
from ThackTech import bdgtools



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('wig', required=True, help='WIG file to convert into bedgraph.')
    parser.add_argument('-o', '--output', help='Specifies the output file. Default is stdout')
    parser.add_argument('-f', '--format', choices=['bdg', 'bw'], default='bdg', help='Output format.')
    parser.add_argument('-g', '--genome', required=True, action='store', help='Genome chromosome sizes. Can specify UCSC genome builds (i.e. hg19, mm9) or the location of a chromosome sizes file (standard UCSC genome sizes format).')
    parser.add_argument('-r', '--repairoverlaps', action='store_true', help='Repair overlapping intervals in the bedgraph, substituting the intersecting regions with a new interval with an average score.')
    parser.add_argument('-m', '--method', choices=bdgtools.score_funcs.keys(), default='mean', help='Method to use for computation of score when merging overlapping intervals.')
    parser.add_argument('-e', '--missingregions', choices=['ignore', 'zero'], action='store', default='ignore', help='How to treat regions with no data.')
    parser.add_argument('-q', '--quiet', action='store_true', help='Be as quiet as possible.')
    args = parser.parse_args()
    
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
        
    if args.quiet:
        sys.stderr = os.devnull
    
    chrom_sizes = bdgtools.get_chrom_sizes(args.genome)
    
    buff = cStringIO.StringIO()
    bdgtools.wig_to_bedgraph(open(args.wig, 'r'), buff)
    buff.seek(0)
    
    
    if args.repairoverlaps:
        intervals = bdgtools.repair_overlapping_segments(buff, chrom_sizes, bdgtools.score_funcs[args.method])
        buff = cStringIO.StringIO()
        bdgtools.write_bedgraph(intervals, buff)
        buff.seek(0)
    
    if not args.missingregions == 'ignore':
        intervals = bdgtools.get_complement_fast(buff, chrom_sizes, args.missingregions)
        buff = cStringIO.StringIO()
        bdgtools.write_bedgraph(intervals, buff)
        buff.seek(0)
        
    if args.format == 'bw':
        sys.stdout.write(bdgtools.convert_bdg_to_bw(buff.getvalue(), chrom_sizes))
    else:
        sys.stdout.write(buff.getvalue())
    return
#end main()


if __name__ == "__main__":
    main()
