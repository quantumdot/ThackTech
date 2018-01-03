import os
import sys
import argparse
from ThackTech.aligntools import extract_barcode
from ThackTech import filetools


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('fastq', nargs='+', help="fastq files to extract barcodes from")
	args = parser.parse_args()
	
	for fq in args.fastq:
		fq_name = filetools.basename_noext(fq, True)
		sys.stdout.write("{}\t{}\n".format(fq_name, "\t".join(extract_barcode(fq, True))))
	
#end main()


if __name__ == "__main__":
	main()