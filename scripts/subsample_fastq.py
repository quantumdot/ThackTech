#!/usr/bin/env python

import os
import sys
import argparse
import random 
import itertools
import HTSeq




def main():
    parser = argparse.ArgumentParser(description="Samples a random subset of sequences from fastq or fasta formatted files. "
                                                +"Both single-end and paired-end data is supported. "
                                                +"Files may be .gz compressed or uncompressed.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    single_end_group = parser.add_argument_group("Single-end Processing options")                                            
    single_end_group.add_argument("--inreads", help="Use for single end reads.")
    single_end_group.add_argument("--outreads", help="Use for single end reads.")
    
    paired_end_group = parser.add_argument_group("Paired-end Processing options") 
    paired_end_group.add_argument("--inread1", help="Read 1 of paired-end reads.")
    paired_end_group.add_argument("--inread2", help="Read 2 of paired-end reads.")
    paired_end_group.add_argument("--outread1", help="Read 1 of paired-end reads.")
    paired_end_group.add_argument("--outread2", help="Read 2 of paired-end reads.")
    
    parser.add_argument("--pe", action="store_true", help="Are the reads paired?")
    parser.add_argument("--fraction", type=float, default=0.1, help="Fraction of reads to return.")
    args = parser.parse_args()
    
    
    if args.pe:
        if None in (args.inread1, args.inread2, args.outread1, args.outread2):
            parser.error("In paired-end mode, you must specify all of --inread1, --inread2, --outread1, and --outread2")
        
        in1 = iter(HTSeq.FastqReader(args.inread1))
        in2 = iter(HTSeq.FastqReader(args.inread2))
        out1 = open(args.outread1, "w")
        out2 = open(args.outread2, "w")
        
        for read1, read2 in itertools.izip(in1, in2):
            if random.random() < args.fraction:
                read1.write_to_fastq_file(out1)
                read2.write_to_fastq_file(out2)
        
        out1.close()
        out2.close()
        
    else:
        if None in (args.inreads, args.outreads):
            parser.error("In paired-end mode, you must specify all of --inreads, and --outreads")
        
        inreads = iter(HTSeq.FastqReader(args.inreads))
        outreads = open(args.outreads, "w")
        
        for read in inreads:
            if random.random() < args.fraction:
                read.write_to_fastq_file(outreads)
        
        outreads.close()
#end main()  
    
if __name__ == "__main__":
    main()
    
