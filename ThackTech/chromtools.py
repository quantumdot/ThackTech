import subprocess
import pybedtools
import sys
import cStringIO
import csv
import os


class ChromosomeSets:
    def __init__(self, all, common, uncommon, use=[]):
        self.all = set(all)
        self.common = set(common)
        self.uncommon = set(uncommon)
        self.use = set(use)
    #end __init__()
#end class ChromosomeSets

def get_bigwig_chroms(bw_file):
    """Produce the set of unique chromosome names in a bigwig file
    
    requires UCSC program bigWigInfo be on the $PATH
    """
    chroms = subprocess.check_output('bigWigInfo -chroms "'+bw_file+'" | grep -Po "^[\s]+\K([\w]+)"', shell=True)
    #chroms = subprocess.check_output(r'bigWigInfo -chroms "%s" | sed -r "s/^[[:space:]]+([[:alnum:]]+)[[:space:]]+[[:digit:]]+[[:space:]]+([[:digit:]]+)/\1\t\2/;tx;d;:x"' % (bw_file,), shell=True)
    #return set([ChromosomeInfo(c.split('\t')) for c in chroms.splitlines()])
    return set(chroms.splitlines())
#end get_bigwig_chroms()

def get_bed_chroms(bed_file):
    """Produce the set of unique chromosome names in a bed file
    """
    chroms = set([])
    bed = pybedtools.BedTool(bed_file)
    for i in bed:
        chroms.add(i.chrom)
    return chroms
#end get_bed_chroms()

def chrom_to_int(chrom):
    """Converts a chromosome name from (string)chr12 -> (int)12
    """
    ord3 = lambda x : '%.3d' % ord(x)
    return int(''.join(map(ord3, chrom))) 
#end chrom_to_int()

def get_common_chroms(beds, sigs):
    """Produce a ChromosomeSets describing the union, intersection, and the complement of the intersection
    
    """
    common_chroms = None
    all_chroms = set()
    for s in sigs:
        chroms = get_bigwig_chroms(s)
        all_chroms = all_chroms | chroms
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
    return ChromosomeSets(all=all_chroms, common=common_chroms, uncommon=all_chroms - common_chroms, use=[])
#end get_common_chroms()




class ChromSizes(dict):
    """Represents the UCSC Chromosome Sizes file format
    
    The class inherits and acts like a normal python dict, where keys are 
    chromosome names and valuse are the size of the associated chromosome.
    """
    def __init__(self, genome, path=None):
        self.genome = genome
        if os.path.isfile(genome):
            self.path = genome
        else:
            self.path = path
        
        if self.path is None:
            self.__load_UCSC()
        else:
            self.__parse_chrom_sizes()
    #end __init__()

    def write(self, outhandle):
        """Writes this ChromSizes to outhandle (file already opened for writing)
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
        """Parses a chromosome size file from the file path provided
        """
        with open(self.path, 'r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for row in reader:
                self[row[0]] = int(row[1])
    #end parse_chrom_sizes()
    
    def __load_UCSC(self):
        """Loads chromosome size info from UCSC MySQL database
        
        REQUIRES MySQL be setup
        """
        assert self.genome is not None, "Genome must be specified!"
        args = ['mysql', 
                '--user=genome', 
                '--host=genome-mysql.cse.ucsc.edu', 
                '-ABN', 
                '-D', str(self.genome), 
                '-e', 'SELECT chrom, size from chromInfo'
        ]
        #sys.stderr.write(" ".join(args)+"\n")
        infile = subprocess.check_output(args, stderr=sys.stderr)
    
        rows = csv.reader(cStringIO.StringIO(infile), delimiter='\t')
        for row in rows:
            self[row[0]] = int(row[1])
    #end fetch_UCSC()
#end class ChromSizes


class ChromosomeInfo(object):
    def __init__(self, chrom, size=None):
        if hasattr(chrom, '__iter__'):
            self.chr = chrom[0]
            self.size = chrom[1]
        else:
            self.chr = chrom
            self.size = size
    #end __init__()
    def __repr__(self):
        return "ChromosomeInfo(%s, %s)" % (self.chr, self.size)
    def __str__(self):
        return self.chr
    def __ne__(self, other):
        return (not self.__eq__(other))
    def __eq__(self, other):
        if isinstance(other, ChromosomeInfo):
            return (self.chr == other.chr)
        else:
            return False
    def __hash__(self):
        return hash(self.__str__())
    
#end class ChromosomeInfo






