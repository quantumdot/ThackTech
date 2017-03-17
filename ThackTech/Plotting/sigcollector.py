import pybedtools
import numpy as np



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
    __allowed_alignments = ['center', 'left', 'right', 'scale']
    
    def __init__(self, bed, align, upstream=0, downstream=0, direction=True, genome=None, white_chroms=None):
        self.bed = bed
        assert align in self.__allowed_alignments, "Align must be one of {"+", ".join(self.__allowed_alignments)+"}"
        self.align = align
        self.up = upstream
        self.down = downstream
        self.direction = direction
        if not isinstance(genome, ChromSizes):
            self.genome = ChromSizes(genome)
        else:
            self.genome = genome
        self.white_chroms
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
        
        if alignment == 'center':
            sys.stderr.write("-> producing center points...\n")
            data = midpoint_generator(bedtool)
        elif alignment == 'left':
            sys.stderr.write("-> producing 5' points...\n")
            data = five_prime_generator(bedtool)
        elif alignment == 'right':
            sys.stderr.write("-> producing 3' points...\n")
            data = three_prime_generator(bedtool)
        elif alignment == 'scale':
            sys.stderr.write("-> producing scaled regions...\n")
            data = scaled_interval_generator(bedtool)
        
        return data
    #end provide_intervals()
    
    def generate_midpoints(self):
        """Generates intervals aligned on the center of bed intervals
                                    
        """
        for interval in bedtool:
            midpoint = interval.start + (interval.stop - interval.start) / 2
            if self.direction and interval.strand == '-':
                start = midpoint - self.downstream
                stop  = midpoint + self.upstream
            else:
                start = midpoint - self.upstream
                stop  = midpoint + self.downstream
                
            start, stop = self.clamp_coordinates(interval.chrom, start, stop)
            yield pybedtools.cbedtools.Interval(interval.chrom, start, stop, strand=interval.strand, name=interval.name, score=interval.score)
    #end generate_midpoints()
    
    def generate_left(self):
        for interval in bedtool:
            if self.direction and interval.strand == '-':
                start = feature.stop - self.downstream
                stop  = feature.stop + self.upstream
            else:
                start = feature.start - self.upstream
                stop  = feature.start + self.downstream
                
            start, stop = self.clamp_coordinates(interval.chrom, start, stop)
            yield pybedtools.cbedtools.Interval(interval.chrom, start, stop, strand=interval.strand, name=interval.name, score=interval.score)
    #end midpoint_generator()
    
    def generate_right(self):
        for interval in bedtool:
            if self.direction and interval.strand == '-':
                start = feature.start - self.downstream
                stop  = feature.start + self.upstream
            else:
                start = feature.stop - self.upstream
                stop  = feature.stop + self.downstream
                
            start, stop = self.clamp_coordinates(interval.chrom, start, stop)
            yield pybedtools.cbedtools.Interval(interval.chrom, start, stop, strand=interval.strand, name=interval.name, score=interval.score)
    #end midpoint_generator()
    
    def generate_scaled(self):
        for interval in bedtool:
            if interval.strand == '-':
                upstream   = ((interval.start - self.downstream),  interval.start)
                gene_body  = (interval.start, interval.stop)
                downstream = (interval.stop, (interval.stop + self.upstream))
            else:
                upstream   = ((interval.start - self.upstream), interval.start)
                gene_body  = (interval.start, interval.stop)
                downstream = (interval.stop, (interval.stop + self.downstream))
            
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


# def midpoint_generator(bedtool, up, down):
#     for i in bedtool:
#         midpoint = i.start + (i.stop - i.start) / 2
#         if i.strand == '-':
#             start = np.clip(midpoint - down, 0, sys.maxint)
#             stop = midpoint + up
#         else:
#             start = np.clip(midpoint - up, 0, sys.maxint)
#             stop = midpoint + down
#         yield pybedtools.cbedtools.Interval(i.chrom, start, stop, strand=i.strand, name=i.name, score=i.score)
# #end midpoint_generator()
# 
# 
# 
# def three_prime_generator(bedtool, up, down):
#     for interval in bedtool:
#         yield pybedtools.featurefuncs.three_prime(interval, upstream=up, downstream=down)#, genome='mm9')
# #end midpoint_generator()
# 
# def scaled_interval_generator(bedtool, up, down):
#     for gene in bedtool:
#         if gene.strand == '-':
#             upstream   = pybedtools.cbedtools.Interval(gene.chrom, np.clip((gene.start - down), 0, sys.maxint), gene.start,         strand='-', name=gene.name, score=gene.score)
#             gene_body  = pybedtools.cbedtools.Interval(gene.chrom, gene.start,             gene.stop,             strand='-', name=gene.name, score=gene.score)
#             downstream = pybedtools.cbedtools.Interval(gene.chrom, gene.stop,             (gene.stop + up),     strand='-', name=gene.name, score=gene.score)
#         else:
#             upstream   = pybedtools.cbedtools.Interval(gene.chrom, np.clip((gene.start - up), 0, sys.maxint),     gene.start,         strand='+', name=gene.name, score=gene.score)
#             gene_body  = pybedtools.cbedtools.Interval(gene.chrom, gene.start,             gene.stop,             strand='+', name=gene.name, score=gene.score)
#             downstream = pybedtools.cbedtools.Interval(gene.chrom, gene.stop,             (gene.stop + down),    strand='+', name=gene.name, score=gene.score)
#         gene_body.name = gene.name
#         yield (upstream, gene_body, downstream)
# #end scaled_interval_generator()
# 
# def expand_bed(bed, up=0, down=0, alignment='center', direction=True, white_chroms=None):
#     
# #end expand_bed()




def get_signal(bedtool, sig_file, inp_bed, input_file, label):
    cache_dir = os.path.abspath(gopts['args'].cachedir)
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

    cache_name = os.path.join(cache_dir, "%s.%s" % (gopts['output_base'], filetools.make_str_filename_safe(label)))
    label_input = label+'_input'
    
    sys.stderr.write("-> input is '%s'...\n" %  (input_file,))
    
    if not gopts['args'].cache or not os.path.exists(cache_name + '.npz'):
        sys.stderr.write("-> Loading signal....\n")
        sig = metaseq.genomic_signal(sig_file, detect_signal_type(sig_file))
        sys.stderr.write("-> Computing signal at intervals....\n")
        sig_array = sig.array(bedtool, bins=gopts['num_bins'], stranded=gopts['args'].dir, method=gopts['args'].collectionmethod, processes=gopts['args'].cpus, zero_inf=False, zero_nan=False)
        
        if input_file is not None:
            sys.stderr.write("-> Loading input signal....\n")
            inp_sig = metaseq.genomic_signal(input_file, detect_signal_type(input_file))
            sys.stderr.write("-> Computing input signal at intervals....\n")
            input_array = inp_sig.array(inp_bed, bins=gopts['num_bins'], stranded=gopts['args'].dir, method=gopts['args'].collectionmethod, processes=gopts['args'].cpus, zero_inf=False, zero_nan=False)
            
        if gopts['args'].cache:
            sys.stderr.write("-> Persisting data to disk...\n")
            cache_data = {label: sig_array}
            if input_file is not None:
                cache_data[label_input] = input_array
            metaseq.persistence.save_features_and_arrays(features=bedtool,
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


