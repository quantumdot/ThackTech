import subprocess
import pybedtools
import sys
import cStringIO
import csv


class ChromosomeSets:
    def __init__(self, all, common, uncommon, use=None):
        self.all = set(all)
        self.common = set(common)
        self.uncommon = set(uncommon)
        self.use = set(use)
    #end __init__()
#end class ChromosomeSets

def get_bigwig_chroms(bw_file):
    chroms = subprocess.check_output('bigWigInfo -chroms "'+bw_file+'" | grep -Po "^[\s]+\K([\w]+)"', shell=True)
    return set(chroms.splitlines())
#end get_bigwig_chroms()

def get_bed_chroms(bed_file):
    chroms = set([])
    bed = pybedtools.BedTool(bed_file)
    for i in bed:
        chroms.add(i.chrom)
    return chroms
#end get_bed_chroms()

def chrom_to_int(chrom):
    ord3 = lambda x : '%.3d' % ord(x)
    return int(''.join(map(ord3, chrom))) 
#end chrom_to_int()

def get_common_chroms(beds, sigs):
    common_chroms = None
    all_chroms = set()
    for s in sigs:
        chroms = get_bigwig_chroms(s)
        all_chroms = all | chroms
        if common_chroms is None:
            common_chroms = chroms
        else:
            common_chroms = chroms & common_chroms
    for b in beds:
        chroms = get_bed_chroms(b)
        all_chroms = all_chroms | chroms
        if common_chroms is None:
            common_chroms = chroms
        else:
            common_chroms = chroms & common_chroms
    return ChromosomeSets(all=all_chroms, common=common_chroms, uncommon=all-common_chroms, use=None)
#end get_common_chroms()




class ChromSizes(dict):
    
    def __init__(self, genome, path=None):
        self.genome = genome
        self.path = path
        if self.path is None:
            self.__load_UCSC()
        else:
            self.__parse_chrom_sizes()
    #end __init__()

    def write(self, outhandle):
        """Writes this ChromSizes to outhandle
        format is:
        {chromosome}[TAB]{Size}[NEWLINE]
        """
        for chrom, size in self.iteritems():
            outhandle.write('%s\t%d\n' % (chrom, size))
        outhandle.flush()
    #end write()
    
    def save(self, filename):
        """Saves this ChromSizes to the file specified
        format is:
        {chromosome}[TAB]{Size}[NEWLINE]
        """
        with open(filename, 'w+') as f:
            self.write(f)
    #end save()
    
    def __parse_chrom_sizes(self):
        with open(self.path, 'r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for row in reader:
                self[row[0]] = int(row[1])
    #end parse_chrom_sizes()
    
    def __load_UCSC(self):
        infile = subprocess.check_output(['mysql', 
                                          '--user=genome', 
                                          '--host=genome-mysql.cse.ucsc.edu', 
                                          '-ABN', 
                                          '-D', self.genome, 
                                          '-e', 'SELECT chrom, size from chromInfo'
                                        ], stderr=sys.stderr)
    
        rows = csv.reader(cStringIO.StringIO(infile), delimiter='\t')
        for row in rows:
            self[row[0]] = int(row[1])
    #end fetch_UCSC()
#end class ChromSizes






