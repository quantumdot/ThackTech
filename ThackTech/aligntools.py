import re
import pysam


def is_bam_PE(bam_path):
    '''
        Find out if a sam or bam file contains even one paired-end read
    '''
    with pysam.AlignmentFile(bam_path, "rb") as samfile:
        for read in samfile.fetch():
            if read.is_paired:
                return True
    return False
#end is_bam_PE()


def check_barcode(fastq_path, barcode):
    ''' Check barcode value in a fastq file    
    '''
    return barcode in extract_barcode(fastq_path)
#end check_barcode()

def extract_barcode(fastq_path, complete=False):
    '''Attempts to extract barcode/index sequence(s) from a FASTQ file
    
        This currently only really works for fastq files generated using casava
        versions 1.4 through 1.8 (or at least versions that append the index sequence
        to the read name)
        
        examples that are detected:
        @HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1
        @HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG
        @SN7001246:318:H2W2MBCX2:1:1105:1992:1992 1:N:0:TAGCTT
        
        Parameters:
            fastq_path: path to fastq file to interrogate
            complete: (bool) if true, return the set of all detected barcodes, otherwise only return the first detected barcode
            
        Returns: 
            (set) Set of detected barcode(s)
    '''
    illumina14_18 = re.compile(r"(?:#|:)([AGCT]+)(?:/?)", re.I)
    all_barcodes = set()
    with pysam.FastxFile(fastq_path) as fh:
        for entry in fh:
            print entry.comment
            match = illumina14_18.search(entry.comment)
            if match:
                print "Match found! {}".format(match.group(1).upper())
                all_barcodes.add(match.group(1).upper())
                if not complete:
                    break    
    return all_barcodes
#end extract_barcode()
                
    