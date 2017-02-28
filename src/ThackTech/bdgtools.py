import re
import sys
import tempfile
import subprocess
import numpy as np
from ThackTech import chromtools
from functools import total_ordering
from intervaltree import IntervalTree
from ThackTech.chromtools import ChromSizes


score_funcs = {}
score_funcs['mean'] 	= np.mean
score_funcs['median'] 	= np.median
score_funcs['sum'] 		= np.sum
score_funcs['min'] 		= np.amin
score_funcs['max'] 		= np.amax
score_funcs['absmin'] 	= lambda lst: min([abs(i) for i in lst])
score_funcs['absmax'] 	= lambda lst: max([abs(i) for i in lst])
score_funcs['count'] 	= lambda lst: len(lst)


@total_ordering
class BedGraphInterval:
	def __init__(self, chrom, start, stop, score):
		self.chr = chrom
		self.start = start
		self.stop = stop
		self.score = score
	#end BedGraphInterval.__init__()
	
	@property
	def numeric_chrom(self):
		return chromtools.chrom_to_int(self.chr)
	#end numeric_chrom()
	
	def getOverlap(self, other):
		if self.chr != other.chr:
			return None
		if self.stop < other.start:
			return None
		if self.stop >= other.start:
			return BedGraphInterval(self.chr, other.start, self.stop, 0)
	#end Interval.getOverlap()
	
	def __radd__(self, other):
		return other + self.score
	
	def __abs__(self):
		return abs(self.score)
	
	def __str__(self):
		return "%s\t%d\t%d\t%f" % (self.chr, self.start, self.stop, self.score)
	#end __str__()
	
	def __repr__(self):
		return "BedGraphInterval(%s)" % (self.chr, self.start, self.stop, self.score)
	#end __repr__()
	
	def __eq__(self, other):
		return ((self.chr.lower(), self.start, self.stop) == (other.chr.lower(), other.start, other.stop))

	def __lt__(self, other):
		return ((self.chr.lower(), self.start, self.stop) < (other.chr.lower(), other.start, other.stop))
#end class Interval

def __parse_part(line, regex, default=None, group=1):
	match = re.search(regex, line)
	if match is not None:
		return match.group(group)
	else:
		return default
#end parse_part()

def parse_bedgraph(input_str, sort=True):
	"""Parses a bedgraph formatted string into a list of BedGraphIntervals
	
	"""
	results = []
	for line in input_str:
		match = re.search(r"(.+)\s+(\d+)\s+(\d+)\s+(-?[\d\.e-]+)", line)
		if match is not None:
			results.append(BedGraphInterval(match.group(1), int(match.group(2)), int(match.group(3)), float(match.group(4))))
		else:
			sys.stderr.write('WARNING: Found non-conforming line: %s\n' % (line,))
	if sort:
		results.sort()
	return results
#end parse_bedgraph()

def parse_wig(input_str, sort=True):
	"""Parses a wig formatted string into a list of BedGraphIntervals
	
	Accepts variable- or fixed-step wig files
	TODO: add UCSC link to wig specification
	"""
	mode = ''
	currChrom = ''
	lastpos = 0
	span = 0
	step = 0
	results = []
	for line in input_str:
		line = line.strip()
		if line.startswith('#') or line.startswith('browser') or line.startswith('track'):
			continue # we ignore these types of lines....
		elif line.startswith('variableStep') or line.startswith('fixedStep'):
			mode = __parse_part(line, r"(variableStep|fixedStep)")
			currChrom = __parse_part(line, r"chrom=([\w\d]+)")
			span = int(__parse_part(line, r"span=(\d+)", 1))
			step = int(__parse_part(line, r"step=(\d+)", 0))
			lastpos = int(__parse_part(line, r"start=(\d+)", 0))
			sys.stderr.write('Processing %s in %s mode\n' % (currChrom, mode))
			continue
		else:
			if mode == 'variableStep':
				parts = re.match(r"^(\d+)\s([\d\.-e]+)", line)
				pos = int(parts.group(1))
				value = parts.group(2)
				results.append(BedGraphInterval(currChrom, pos-1, pos+span, value))
				
			elif mode == 'fixedStep':
				parts = re.match(r"^([\d\.-e]+)", line)
				value = parts.group(1)
				results.append(BedGraphInterval(currChrom, lastpos-1, lastpos+span, value))
				lastpos += step
	if sort:
		results.sort()
	return results
#end wig_to_bedgraph()

def parse_bed(input_str, sort=True):
	"""Parses a bed formatted string into a list of BedGraphIntervals
	
	Accepts score column is interpreted as BedGraphInterval score. Strand information is ignored.
	TODO: add UCSC link to bed specification
	"""
	results = []
	for line in input_str:
		line = line.strip()
		if line.startswith('#') or line.startswith('browser') or line.startswith('track'):
			continue # we ignore these types of lines....
		else:
			parts = re.match(r"^(\w+)\s+(\d+)\s+(\d+)\s+(?:.+?)\s+(-?[\.\de-]+)", line)
			chrom = parts.group(1)
			start = int(parts.group(2))
			stop = int(parts.group(3))
			value = parts.group(4)
			results.append(BedGraphInterval(chrom, start, stop, value))
	if sort:
		results.sort()
	return results
#end wig_to_bedgraph()

def write_bedgraph(intervals, outhandle):
	"""Writes a list of BedGraphIntervals to outhandle in bedgraph format
	
	Parameters:
		intervals: (list of BedGraphIntervals) Intervals to write as bedgraph
		outhandle: (file handle) file handle to write to
		
	"""
	for interval in intervals:
		outhandle.write("%s\t%d\t%d\t%f\n" % (interval.chr, interval.start, interval.stop, interval.score))
	outhandle.flush()
#end write_bedgraph()

def write_bigwig(intervals, chrsizes, outhandle):
	"""Writes a list of BedGraphIntervals to outhandle in BigWig format
	
	Parameters:
		intervals: (list of BedGraphIntervals) Intervals to write to BigWig
		chrsizes: (chromtools.ChromSizes) Chromosome size information
		outhandle: (file handle) file handle to write to
		
	This method depends on the program `bedGraphToBigWig` being on the PATH
	"""
	with tempfile.NamedTemporaryFile() as temp_bdg:
		write_bedgraph(intervals, temp_bdg)
		with tempfile.NamedTemporaryFile() as temp_chrsizes:
			chrsizes.write(temp_chrsizes)
			temp_chrsizes.flush()
			with tempfile.NamedTemporaryFile() as temp_bw:
				p = subprocess.Popen(['bedGraphToBigWig', temp_bdg.name, temp_chrsizes.name, temp_bw.name], stderr=sys.stderr)
				p.communicate()
				temp_bw.seek(0)
				outhandle.write(temp_bw.read())
	outhandle.flush()
#end convert_bdg_to_bw()

def clip_chrom_sizes(intervals, chromsizes):
	"""Clip intervals to be within the size bounds defined by chromsizes
	
	TODO: complete implementation!
	
	Parameters:
		intervals: (list of BedGraphInterval) intervals to repair
		chrsizes: (chromtools.ChromSizes) Chromosome size information
		
	Returns:
		list of BedGraphIntervals conforming to chromosome sizes
	"""
	trees = {}
	#init a intervaltree for each chromosome
	for chrom in chromsizes:
		trees[chrom] = IntervalTree()
	
	#populate the chromosome-specific intervaltrees with the intervals
	for i in range(len(intervals)):
		trees[intervals[i].chr].addi(intervals[i].start, intervals[i].stop, intervals[i].score)
		
		
#end clip_chrom_sizes()


def reduce_overlaps(intervals, chromsizes, scorefunc):
	"""Resolves overlapping BedGraphIntervals by reducing overlapping values by `scorefunc`
	
	Parameters:
		intervals: (list of BedGraphInterval) intervals to repair
		chrsizes: (chromtools.ChromSizes) Chromosome size information
		scorefunc: (callable) Reduce function, passed list of overlapping values, return new value
	
	Returns:
		list of BedGraphIntervals with overlaps resolved
	"""
	sys.stderr.write('Repairing overlapping regions\n')
	
	trees = {}
	#init a intervaltree for each chromosome
	for chrom in chromsizes:
		trees[chrom] = IntervalTree()
	
	#populate the chromosome-specific intervaltrees with the intervals
	for i in range(len(intervals)):
		trees[intervals[i].chr].addi(intervals[i].start, intervals[i].stop, intervals[i].score)

	#split overlaps, and then merge the the now equal-coordinate intervals,
	#accumulating their values along the way as a list
	results = []
	for chrom in chromsizes:
		sys.stderr.write(' -> working on %s....\n' % (chrom,))
		trees[chrom].split_overlaps()
		trees[chrom].merge_equals(lambda l, v: l + [v], [])
	
		for i in trees[chrom]:
			results.append(BedGraphInterval(chrom, i.begin, i.end, scorefunc(i.data)))
			
	results.sort()
	return results
#end reduce_overlaps()

def fill_complement(intervals, chroms, fill_value):
	"""fills the complement of the bedgraph with some value
	
	Args:
		f (string): b
	"""
	sys.stderr.write('Finding regions with missing data\n')
	
	genome_trees = {}
	for chrom in chroms:
		genome_trees[chrom] = IntervalTree()
		genome_trees[chrom].addi(0, chroms[chrom], fill_value)

	
	results = []
	for iv in intervals:
		results.append(iv)
		genome_trees[iv.chr].chop(iv.start, iv.stop)
		
	for chrom in chroms:
		for iv in genome_trees[chrom]:
			results.append(BedGraphInterval(chrom, iv.begin, iv.end, iv.data))	
			
	results.sort()
	return results
#end fill_complement()

