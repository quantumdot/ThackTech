import os
import sys
import pybedtools
from pybedtools import Interval
import numpy as np
import metaseq
import subprocess
from ThackTech import filetools
from ThackTech import chromtools



def detect_signal_type(sig_file):
    """Given a filename, attempts to determine the file format based on extension
    
    """
    ext = os.path.splitext(sig_file)[1][1:].strip().lower()
    if ext in ['bigwig', 'bw']:
        return "bigwig"
    else:
        return ext
#end detect_signal_type()


class IntervalProvider:
    
    
    def __init__(self, bedfile, collector_opts, genome=None, white_chroms=None):
        self.bed = bedfile
        self.co = collector_opts
        self.white_chroms = white_chroms
        
        if not isinstance(genome, chromtools.ChromSizes):
            self.genome = chromtools.ChromSizes(genome)
        else:
            self.genome = genome
        
    #end __init__()
    
    def clamp_coordinates(self, chrom, start, stop):
        return (start, stop)
    
    
        if self.genome is not None:
            stop = min(stop, self.genome[chrom])
            start = max(start, 0)
        if start < 0:
            start = 0
        if start > stop:
            start = stop
        return (start, stop)
    #end clamp_coordinates()
    
    def provide_intervals(self):
        data = None
        
        if self.co.align == 'center':
            sys.stderr.write("-> producing center points...\n")
            data = self.generate_midpoints(self.generate_origional())
        elif self.co.align == 'left':
            sys.stderr.write("-> producing 5' points...\n")
            data = self.generate_left(self.generate_origional())
        elif self.co.align == 'right':
            sys.stderr.write("-> producing 3' points...\n")
            data = self.generate_right(self.generate_origional())
        elif self.co.align == 'scale':
            sys.stderr.write("-> producing scaled regions...\n")
            data = self.generate_scaled(self.generate_origional())
        elif self.co.align == 'origional':
            sys.stderr.write("-> producing origional regions...\n")
            data = self.generate_origional()
        
        return data
    #end provide_intervals()
    
    def generate_origional(self):
        """Generates the origional intervals this provider wraps
        """
        bedtool = pybedtools.BedTool(self.bed)
        if self.white_chroms is not None:
            bedtool = bedtool.filter(lambda d: d.chrom in self.white_chroms)
        
        for interval in bedtool:
            yield interval
    #end generate_midpoints()
    
    def generate_origional_as_tuples(self):
        """Generates the origional intervals this provider wraps
        """
        bedtool = pybedtools.BedTool(self.bed)
        if self.white_chroms is not None:
            bedtool = bedtool.filter(lambda d: d.chrom in self.white_chroms)
        
        for i in bedtool:
            yield (i.chrom, i.start, i.stop, i.name, i.score, i.strand)
    #end generate_midpoints()
    
    def generate_midpoints(self, bedtool):
        """Generates intervals aligned on the center of bed intervals
                                    
        """
        for interval in bedtool:
            midpoint = interval.start + (interval.stop - interval.start) / 2
            if self.co.direction and interval.strand == '-':
                start = midpoint - self.co.downstream
                stop  = midpoint + self.co.upstream
            else:
                start = midpoint - self.co.upstream
                stop  = midpoint + self.co.downstream
                
            start, stop = self.clamp_coordinates(interval.chrom, start, stop)
            yield Interval(interval.chrom, start, stop, strand=interval.strand, name=interval.name, score=interval.score)            
    #end generate_midpoints()
    
    def generate_left(self, bedtool):
        for interval in bedtool:
            if self.co.direction and interval.strand == '-':
                start = interval.stop - self.co.downstream
                stop  = interval.stop + self.co.upstream
            else:
                start = interval.start - self.co.upstream
                stop  = interval.start + self.co.downstream
                
            start, stop = self.clamp_coordinates(interval.chrom, start, stop)
            yield Interval(interval.chrom, start, stop, strand=interval.strand, name=interval.name, score=interval.score)
    #end midpoint_generator()
    
    def generate_right(self, bedtool):
        for interval in bedtool:
            if self.co.direction and interval.strand == '-':
                start = interval.start - self.co.downstream
                stop  = interval.start + self.co.upstream
            else:
                start = interval.stop - self.co.upstream
                stop  = interval.stop + self.co.downstream
                
            start, stop = self.clamp_coordinates(interval.chrom, start, stop)
            yield Interval(interval.chrom, start, stop, strand=interval.strand, name=interval.name, score=interval.score)
    #end midpoint_generator()
    
    def generate_scaled(self, bedtool):
        for interval in bedtool:
            if self.co.direction and interval.strand == '-':
                upstream   = ((interval.start - self.co.downstream),  interval.start)
                gene_body  = (interval.start, interval.stop)
                downstream = (interval.stop, (interval.stop + self.co.upstream))
            else:
                upstream   = ((interval.start - self.co.upstream), interval.start)
                gene_body  = (interval.start, interval.stop)
                downstream = (interval.stop, (interval.stop + self.co.downstream))
            
            #make sure the coordinates are sane
            upstream = self.clamp_coordinates(interval.chrom, upstream[0], upstream[1])
            gene_body = self.clamp_coordinates(interval.chrom, gene_body[0], gene_body[1])
            downstream = self.clamp_coordinates(interval.chrom, downstream[0], downstream[1])
            
            #make the actual interval
            upstream   = Interval(interval.chrom, upstream[0], upstream[1], strand=interval.strand, name=interval.name, score=interval.score)
            gene_body  = Interval(interval.chrom, gene_body[0], gene_body[1], strand=interval.strand, name=interval.name, score=interval.score)
            downstream = Interval(interval.chrom, downstream[0], downstream[1], strand=interval.strand, name=interval.name, score=interval.score)
            
            #pybedtools screws up when creating a standalone interval and does not set the name correctly from the constructor
            #re-set the name manually!
            upstream.name = interval.name
            gene_body.name = interval.name
            downstream.name = interval.name
            yield (upstream, gene_body, downstream)
    #end generate_scaled()
#end class IntervalProvider




class CollectorOptions:
    
    def __init__(self):
        self.align = 'center'
        self.upstream = 0
        self.downstream = 0
        self.resolution = 1
        self.direction = True
        self.scaleregionsize = 3000
    #end __init()
    
    __allowed_alignments = ['center', 'left', 'right', 'scale']
    def validate(self):
        """Validate this CollectorOptions to ensure validity
        """
        assert self.align in self.__allowed_alignments, "Align must be one of {"+", ".join(self.__allowed_alignments)+"}"
        assert self.upstream >= 0, "Upstream must be >= 0!"
        assert self.downstream >= 0, "Downstream must be >= 0!"
        assert self.resolution > 0, "Resolution must be > 0!"
        assert isinstance(self.direction, bool), "Direction must be boolean!"
        assert self.scaleregionsize > 0, "Scale region size must be >0!"
    #end validate()
    
    @property
    def units(self):
        """Get a tuple that describes the units for this collector
        
        The returned tuple contains the following information:
            [0] => denomenator to use for division of the region
            [1] => human readable string describing the units
        
        example: (1e3, 'Kb')
        """
        if (self.total_bins * self.resolution) >= 2e9:
            return (1e9, 'Gb')
        if (self.total_bins * self.resolution) >= 2e6:
            return (1e6, 'Mb')
        elif (self.total_bins * self.resolution) >= 2e3:
            return (1e3, 'Kb')
        else:
            return (1, 'bp')
    #end units()
    
    @property
    def xaxis(self):
        """Retrun a numpy array appropriate for use as an x-axis for plotting based on this CollectorOptions
        """
        if self.align == 'scale':
            return np.linspace((-self.upstream), abs(self.scaleregionsize) + abs(self.downstream), self.total_bins)
        else:
            return np.linspace((-self.upstream), self.downstream, self.total_bins)
    #end xaxis()
    
    @property
    def xaxis_label(self):
        """Get a string that is an appropriate x-axis label
        """
        unit_txt = self.units[1]
        if self.align == 'center':
            return "Distance from Center of Element (%s)" % (unit_txt,)
        elif self.align == 'left':
            return "Distance from 5' end of Element (%s)" % (unit_txt,)
        elif self.align == 'right':
            return "Distance from 3' end of Element (%s)" % (unit_txt,)
        elif self.align == 'scale':
            return "Upstream (%s); %d%s of Meta-Element; Downstream (%s)" % (unit_txt, (abs(self.scaleregionsize) / self.units[0]), unit_txt, unit_txt)
    #end xaxis_label
    
    @property
    def total_bins(self):
        """Returns the absolute total number of bins
        
        """
        if self.align == 'scale':
            return sum(self.num_bins)
        else:
            return self.num_bins
    #end total_bins
    
    @property
    def num_bins(self):
        """Returns the number of bins, which may be an int or a tuple of ints
        
        In the case of scaled region mode, a tuple of three ints is retruned representing (left, scaled, right), otherwise an int is returned
        """
        if self.align == 'scale':
            return ((abs(self.upstream) / self.resolution), (abs(self.scaleregionsize) / self.resolution), (abs(self.downstream) / self.resolution))
        else:
            return (abs(self.upstream) + abs(self.downstream)) / self.resolution
    #end num_bins()
    
    def __str__(self):
        """Gets a string representation of this collector options
        """
        return "%du_%dd_%dr_%s%s" % (self.upstream, self.downstream, self.resolution, (self.align+str(self.scaleregionsize) if self.align == 'scale' else self.align), ('_dir' if self.direction else ''))
#end class CollectorOptions



def get_signal(regions, label, sig_file, input_file=None, cache_dir=None, cache_base="", collectionmethod="get_as_array", norm_method='log2', norm_pseudocount=1.0, cpus=1):
    """ Gets a numpy nd-array representing the signal in sig_file at regions
    
        Parameters:
            regions: an IntervalProvider object representing the regions to profile
            label: (string) Name for this signal set
            sig_file: (string) path to a file accepted by metaseq.genomic_signal
            input_file: (string)[optional] path to a input signal for normalization
            cache_dir: (string) to enable cacheing of results, specify a directory for cache files to be placed
            cache_base: (string) basename for cache files, label and collector options are automatically appended
            collectionmethod: (string) signal collection method to use. see metaseq.genomic_signal.array
            cpus: (int) number of parallel processes to use for signal collection
            norm_method: Normalization method to use when input is included. default is log2, choices are {'ratio', 'log2', 'reciprocal_ratio', 'subtract', 'add', 'mean'}
            norm_pseudocount: (float) small number to add to both operands to avoid division by zero errors
            
        Returns:
            numpy nd-array of signal at regions
    
    """
    if cache_dir is not None:
        cache_dir = os.path.abspath(cache_dir)
        filetools.ensure_dir(cache_dir)
        cache_name = os.path.join(cache_dir, "%s.%s.%s" % (cache_base, regions.co, filetools.make_str_filename_safe(label)))
        label_input = label+'_input'
    
    sys.stderr.write("-> input is '%s'...\n" %  (input_file,))
    #pybedtools.BedTool(regions.provide_intervals()).saveas('bedtool_gen.bed')
    if (cache_dir is None) or (not os.path.exists(cache_name + '.npz')):
        sys.stderr.write("-> Loading signal....\n")
        sig = metaseq.genomic_signal(sig_file, detect_signal_type(sig_file))
        sys.stderr.write("-> Computing signal at intervals....\n")
        #print list(regions.provide_intervals())
        #print regions.co
        #print regions.co.num_bins
        #print list(regions.provide_intervals())[0:5]
        
        sig_array = sig.array(regions.provide_intervals(), bins=regions.co.num_bins, stranded=regions.co.direction, method=collectionmethod, processes=cpus, zero_inf=False, zero_nan=False)
        
        if input_file is not None:
            sys.stderr.write("-> Loading input signal....\n")
            inp_sig = metaseq.genomic_signal(input_file, detect_signal_type(input_file))
            sys.stderr.write("-> Computing input signal at intervals....\n")
            input_array = inp_sig.array(regions.provide_intervals(), bins=regions.co.num_bins, stranded=regions.co.direction, method=collectionmethod, processes=cpus, zero_inf=False, zero_nan=False)
            
        if cache_dir is not None:
            sys.stderr.write("-> Persisting data to disk...\n")
            cache_data = {label: sig_array}
            if input_file is not None:
                cache_data[label_input] = input_array
            metaseq.persistence.save_features_and_arrays(features=[],#regions.provide_intervals(),
                                                         arrays=cache_data,
                                                         prefix=cache_name,
                                                         #link_features=True,
                                                         overwrite=True)
    else:
        sys.stderr.write("-> Loding data from cache....\n")
        features, arrays = metaseq.persistence.load_features_and_arrays(prefix=cache_name)
        sig_array = arrays[label]
        if input_file is not None:
            input_array = arrays[label_input]
            
    if input_file is not None and norm_method is not None:    
        sys.stderr.write("-> Normalizing signal to input....\n")
        if norm_method in ['ratio', 'log2', 'reciprocal_ratio']:
            sig_array = (sig_array + norm_pseudocount) / (input_array + norm_pseudocount)
            if norm_method == 'log2':
                sig_array = np.log2(sig_array)
            elif norm_method == 'reciprocal_ratio':
                # the reciprocal ratio of a/b
                # is a/b if a/b > 1 else -1* b/a
                sig_array = np.divide(-1.0, sig_array, where=(sig_array < 1))
        else:
            if norm_method == 'subtract':
                sig_array = sig_array - input_array
            elif norm_method == 'add':
                sig_array = sig_array + input_array
            elif norm_method == 'mean':
                sig_array = (sig_array + input_array) / 2.0      
    
    return sig_array
#end get_signal()


def get_bed_score_signal(regions):
    matrix = []
    matrix_size = regions.co.total_bins
    for el in regions.generate_origional():
        matrix.append([float(el.score)] * matrix_size)
    return np.array(matrix)
#end get_bed_score_signal()



def get_bed_score_signal_complex(regions, cache_dir=None, cache_base="", collectionmethod="get_as_array", cpus=1, **kwargs):
    if cache_dir is not None:
        cache_dir = os.path.abspath(cache_dir)
        filetools.ensure_dir(cache_dir)
    
    bed_basename = os.path.splitext(os.path.basename(regions.bed))[0]
    bw_name = os.path.join(cache_dir, bed_basename+'.bw')
    
    if not os.path.exists(bw_name):
        cmd = [
            'cleanBedGraph.py',
            '--informat', 'bed', 
            #'--quiet', 
            '--output', bw_name,
            '--outformat', 'bw', 
            '--genome', regions.genome.genome,
            '--clipsize',
            '--repairoverlaps', 'mean', 
            '--missing', '0', 
            regions.bed
        ]
        p = subprocess.Popen(cmd)
        p.communicate()
  
    return get_signal(regions, bed_basename+'_signal', bw_name, None, cache_dir, cache_base, collectionmethod, cpus)
#end get_bed_score_signal()



