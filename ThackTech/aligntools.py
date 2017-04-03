import pysam


def is_bam_PE(bam_path):
    '''
        Find is a sam or bam file contains even one paired-end read
    '''
    with pysam.AlignmentFile(bam_path, "rb") as samfile:
        for read in samfile.fetch():
            if read.is_paired:
                return True
    return False
#end is_bam_PE()