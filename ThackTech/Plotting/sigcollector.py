import pybedtools
import numpy as np
import metaseq
from ThackTech import filetools



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
        self.bed = bed
        self.co = collector_opts
        self.white_chroms
        
        if not isinstance(genome, ChromSizes):
            self.genome = ChromSizes(genome)
        else:
            self.genome = genome
        
    #end __init__()
    
    def clamp_coordinates(self, chrom, start, stop):
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
        bedtool = pybedtools.BedTool(self.bed)
        if white_chroms is not None:
            bedtool = bedtool.filter(lambda d: d.chrom in white_chroms)
        
        if self.co.align == 'center':
            sys.stderr.write("-> producing center points...\n")
            data = self.generate_midpoints(bedtool)
        elif self.co.align == 'left':
            sys.stderr.write("-> producing 5' points...\n")
            data = self.generate_left(bedtool)
        elif self.co.align == 'right':
            sys.stderr.write("-> producing 3' points...\n")
            data = self.generate_right(bedtool)
        elif self.co.align == 'scale':
            sys.stderr.write("-> producing scaled regions...\n")
            data = self.generate_scaled(bedtool)
        elif self.co.align == 'origional':
            sys.stderr.write("-> producing origional regions...\n")
            data = self.generate_origional(bedtool)
        
        return data
    #end provide_intervals()
    
    def generate_origional(self, bedtool):
        for interval in bedtool:
            yield interval
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
            yield pybedtools.cbedtools.Interval(interval.chrom, start, stop, strand=interval.strand, name=interval.name, score=interval.score)
    #end generate_midpoints()
    
    def generate_left(self, bedtool):
        for interval in bedtool:
            if self.co.direction and interval.strand == '-':
                start = feature.stop - self.co.downstream
                stop  = feature.stop + self.co.upstream
            else:
                start = feature.start - self.co.upstream
                stop  = feature.start + self.co.downstream
                
            start, stop = self.clamp_coordinates(interval.chrom, start, stop)
            yield pybedtools.cbedtools.Interval(interval.chrom, start, stop, strand=interval.strand, name=interval.name, score=interval.score)
    #end midpoint_generator()
    
    def generate_right(self, bedtool):
        for interval in bedtool:
            if self.co.direction and interval.strand == '-':
                start = feature.start - self.co.downstream
                stop  = feature.start + self.co.upstream
            else:
                start = feature.stop - self.co.upstream
                stop  = feature.stop + self.co.downstream
                
            start, stop = self.clamp_coordinates(interval.chrom, start, stop)
            yield pybedtools.cbedtools.Interval(interval.chrom, start, stop, strand=interval.strand, name=interval.name, score=interval.score)
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
            upstream   = pybedtools.cbedtools.Interval(interval.chrom, upstream[0], upstream[1], strand=interval.strand, name=gene.name, score=gene.score)
            gene_body  = pybedtools.cbedtools.Interval(interval.chrom, gene_body[0], gene_body[1], strand=interval.strand, name=gene.name, score=gene.score)
            downstream = pybedtools.cbedtools.Interval(interval.chrom, downstream[0], downstream[1], strand=interval.strand, name=gene.name, score=gene.score)
            
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
        if (self.total_bins * self.resolution) > 2e9:
            return (1e9, 'Gb')
        if (self.total_bins * self.resolution) > 2e6:
            return (1e6, 'Mb')
        elif (self.total_bins * self.resolution) > 2e3:
            return (1e3, 'Kb')
        else:
            return (1, 'bp')
    #end units()
    
    @property
    def xaxis(self):
        if self.align == 'scale':
            return np.linspace((-self.upstream), abs(self.scaleregionsize) + abs(self.downstream), self.total_bins)
        else:
            return np.linspace((-self.upstream), self.downstream, self.total_bins)
    #end xaxis()
    
    @property
    def xaxis_label(self):
        unit_txt = self.units[1]
        if args.align == 'center':
            return "Distance from Center of Element (%s)" % (unit_txt,)
        elif args.align == 'left':
            return "Distance from 5' end of Element (%s)" % (unit_txt,)
        elif args.align == 'right':
            return "Distance from 3' end of Element (%s)" % (unit_txt,)
        elif args.align == 'scale':
            return "Upstream (%s); %d%s of Meta-Element; Downstream (%s)" % (unit_txt, (abs(self.scaleregionsize) / self.units[0]), unit_txt, unit_txt)
    #end xaxis_label
    
    @property
    def total_bins(self):
        """Retruns the absolute total number of bins
        
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
            return (abs(self.upstream) / self.resolution, abs(self.scaleregionsize) / self.resolution, abs(self.downstream) / self.resolution)
        else:
            return (abs(self.upstream) + abs(self.downstream)) / self.resolution
    #end num_bins()
    
    def __str__(self):
        """Gets a string representation of this collector options
        """
        return "%du_%dd_%dr_%s%s" % (self.upstream, self.downstream, self.resolution, (self.align+str(self.scaleregionsize) if self.align == 'scale' else self.align), ('_dir' if self.direction else ''))
#end class CollectorOptions



def get_signal(regions, label, sig_file, input_file=None, cache_dir=None, cache_base=""):
    if cache_dir is not None:
        cache_dir = os.path.abspath(cache_dir)
        filetools.ensure_dir(cache_dir)
        cache_name = os.path.join(cache_dir, "%s.%s.%s" % (cache_base, regions.co, filetools.make_str_filename_safe(label)))
        label_input = label+'_input'
    
    sys.stderr.write("-> input is '%s'...\n" %  (input_file,))
    
    if (cache_dir is None) or (not os.path.exists(cache_name + '.npz')):
        sys.stderr.write("-> Loading signal....\n")
        sig = metaseq.genomic_signal(sig_file, detect_signal_type(sig_file))
        sys.stderr.write("-> Computing signal at intervals....\n")
        sig_array = sig.array(regions.provide_intervals(), bins=regions.co.num_bins, stranded=regions.co.direction, method=gopts['args'].collectionmethod, processes=gopts['args'].cpus, zero_inf=False, zero_nan=False)
        
        if input_file is not None:
            sys.stderr.write("-> Loading input signal....\n")
            inp_sig = metaseq.genomic_signal(input_file, detect_signal_type(input_file))
            sys.stderr.write("-> Computing input signal at intervals....\n")
            input_array = inp_sig.array(regions.provide_intervals(), bins=regions.co.num_bins, stranded=regions.co.direction, method=gopts['args'].collectionmethod, processes=gopts['args'].cpus, zero_inf=False, zero_nan=False)
            
        if cache_dir is not None:
            sys.stderr.write("-> Persisting data to disk...\n")
            cache_data = {label: sig_array}
            if input_file is not None:
                cache_data[label_input] = input_array
            metaseq.persistence.save_features_and_arrays(features=regions.provide_intervals(),
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
            
    if input_file is not None:    
        sys.stderr.write("-> Normalizing signal to input....\n")
        sig_array = sig_array - input_array    
    
    return sig_array
#end get_signal()


def get_bed_score_signal(bed, white_chroms=None):
    elements = pybedtools.BedTool(bed)
    if white_chroms is not None:
        elements = elements.filter(lambda d: d.chrom in white_chroms)
    matrix = []
    for el in elements:
        matrix.append([float(el.score)]*gopts['total_bins'])
    return np.array(matrix)
#end get_bed_score_signal()

def get_bed_score_signal_complex(bed, genome, white_chroms=None):
    import subprocess
    cache_dir = os.path.abspath(gopts['args'].cachedir)
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)
    bed_basename = os.path.splitext(os.path.basename(bed))[0]
    bw_name = os.path.join(cache_dir, bed_basename+'.bw')
    if not os.path.exists(bw_name):
        cmd = ['python', '/home/josh/scripts/bedToBedGraph.py', '--quiet', '--output', bw_name, '--genome', genome, '--format', 'bw', '--repairoverlaps', '--method', 'mean', '--missingregions', 'zero', bed]
        #print " ".join(cmd)
        p = subprocess.Popen(cmd)
        p.communicate()

    bedtool1 = expand_bed(gopts['args'].up, gopts['args'].down, gopts['args'].align, bed, gopts['chromsets'].use)        
    signal = get_signal(bedtool1, bw_name, None, None, bed_basename+'_signal')
    return signal
#end get_bed_score_signal()


