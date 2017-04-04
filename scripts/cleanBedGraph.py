#!/usr/bin/env python
import os
import sys
import argparse
from ThackTech import bdgtools, chromtools



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('interval-data', help='File containing interval-based score data to clean.')
    parser.add_argument('-if', '--informat', default='auto', choices=['bdg', 'bed', 'wig', 'bw', 'auto'], help="Input format is auto-detected based on file extension, but use this option to force.")
    parser.add_argument('-o', '--output', help='Specifies the output file. Default is stdout')
    parser.add_argument('-of', '--outformat', choices=['bdg', 'bw'], default='bdg', help='Output format.')
    parser.add_argument('-g', '--genome', action='store', required=True, help='Genome chromosome sizes. Can specify UCSC genome builds (i.e. hg19, mm9) or the location of a chromosome sizes file (standard UCSC genome sizes format).')
    parser.add_argument('-c', '--clipsize', action='store_true', help="Clips the interval data so that it conforms to the chrom sizes of --genome.")
    parser.add_argument('-r', '--repairoverlaps', default=None, help='Repair overlapping intervals in the bedgraph, substituting the intersecting regions with a new interval having a score computed by one of '+str(bdgtools.score_funcs.keys()))
    parser.add_argument('-m', '--missing', default=None, help='Set regions with no data to have some value. Accepted values are a literal number (i.e. 0 or 1.5) or one of '+str(bdgtools.score_funcs.keys())+' to dynamically compute new score from adjacent intervals.')
    parser.add_argument('-q', '--quiet', action='store_true', help='Be as quiet as possible.')
    args = parser.parse_args()
    
    if args.repairoverlaps is not None and args.repairoverlaps not in bdgtools.score_funcs.keys():
        parser.error("--repairoverlaps must be one of {}".format(bdgtools.score_funcs.keys()))
    
        
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
        
    if args.quiet:
        sys.stderr = os.devnull
        
    chrom_sizes = chromtools.ChromSizes(args.genome)
    intervals = bdgtools.open_file_as_bedgraph(args.interval_data, args.informat)
    
    
    if args.repairoverlaps is not None:
        intervals = bdgtools.reduce_overlaps(intervals, chrom_sizes, bdgtools.score_funcs[args.repairoverlaps])
    
    
    if args.missing is not None:
        if args.missing in bdgtools.score_funcs.keys():
            missing_score = bdgtools.score_funcs[args.missing]
        else:
            missing_score = lambda lst: float(args.missing)
        intervals = bdgtools.fill_complement(intervals, chrom_sizes, missing_score)
        

    if args.clipsize:
        intervals = bdgtools.clip_chrom_sizes(intervals, chrom_sizes)
    
    
    if args.outformat == 'bw':
        bdgtools.write_bigwig(intervals, chrom_sizes, sys.stdout)
    else:
        bdgtools.write_bedgraph(intervals, sys.stdout)
    return
#end main()


if __name__ == "__main__":
    main()
