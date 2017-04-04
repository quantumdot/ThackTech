#!/usr/bin/env python

import argparse
import random 
import itertools
import gzip
import pysam



def main():
    parser = argparse.ArgumentParser(description="Samples a random subset of sequences from fastq or fasta formatted files. "
                                                +"Both single-end and paired-end data is supported. Running time is ~O(N) where N=size of file (# of reads). "
                                                +"Files may be .gz compressed or uncompressed."
                                                +"If outfiles have .gz extension data will be gzipped.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    
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
        
        in1 = pysam.FastxFile(args.inread1)
        in2 = pysam.FastxFile(args.inread2)
        
        out1 = open_file_write(args.outread1)
        out2 = open_file_write(args.outread2)
        
        for read1, read2 in itertools.izip(in1, in2):
            if random.random() < args.fraction:
                write_fastq(out1, read1)
                write_fastq(out2, read2)
        
        out1.close()
        out2.close()
        
    else:
        if None in (args.inreads, args.outreads):
            parser.error("In paired-end mode, you must specify all of --inreads, and --outreads")
        
        inreads = pysam.FastxFile(args.inreads)
        outreads = open_file_write(args.outreads)
        
        for read in inreads:
            if random.random() < args.fraction:
                write_fastq(outreads, read)
        
        outreads.close()
#end main()

def open_file_write(filename):
    if filename.lower().endswith((".gz", ".gzip")):
        return gzip.open(filename, "wb")
    else:
        return open(filename, "w")
#end open_file()

def write_fastq(f, read):
    fastq_format_str = "@{name}\n{sequence}\n+\n{quality}\n"
    f.write(fastq_format_str.format(name=read.name, sequence=read.sequence, quality=read.quality))
#end write_fastq()
    
if __name__ == "__main__":
    main()
    
