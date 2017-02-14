#!/usr/bin/python
import re
import sys
import cStringIO
import csv
import numpy as np
import os.path
from ThackTech import chromtools
from intervaltree import IntervalTree
import subprocess
import tempfile

score_funcs = {}
score_funcs['mean'] 	= np.mean
score_funcs['median'] 	= np.median
score_funcs['sum'] 		= np.sum
score_funcs['min'] 		= np.amin
score_funcs['max'] 		= np.amax
score_funcs['absmin'] 	= lambda lst: min([abs(i) for i in lst])
score_funcs['absmax'] 	= lambda lst: max([abs(i) for i in lst])
score_funcs['count'] 	= lambda lst: len(list)




def parse_part(line, regex, default=None, group=1):
	match = re.search(regex, line)
	if match is not None:
		return match.group(group)
	else:
		return default
#end parse_part()


def wig_to_bedgraph(input, output):
	mode = ''
	currChrom = ''
	lastpos = 0
	span = 0
	step = 0
	for line in input:
		if line.strip().startswith('#') or line.startswith('browser') or line.startswith('track'):
			continue # we ignore these types of lines....
		elif line.startswith('variableStep') or line.startswith('fixedStep'):
			mode = parse_part(line, r"(variableStep|fixedStep)")
			currChrom = parse_part(line, r"chrom=([\w\d]+)")
			span = int(parse_part(line, r"span=(\d+)", 1))
			step = int(parse_part(line, r"step=(\d+)", 0))
			lastpos = int(parse_part(line, r"start=(\d+)", 0))
			sys.stderr.write('Processing %s in %s mode\n' % (currChrom, mode))
		else:
			if mode == 'variableStep':
				parts = re.match(r"^(\d+)\s([\d\.-e]+)", line)
				pos = int(parts.group(1))
				value = parts.group(2)
				output.write("%s\t%d\t%d\t%s\n" % (currChrom, pos-1, pos+span, value))
				
			elif mode == 'fixedStep':
				parts = re.match(r"^([\d\.-e]+)", line)
				value = parts.group(1)
				output.write("%s\t%d\t%d\t%s\n" % (currChrom, lastpos-1, lastpos+span, value))
				lastpos += step
#end wig_to_bedgraph()

def bed_to_bedgraph(input, output):
	for line in input:
		if line.strip().startswith('#') or line.startswith('browser') or line.startswith('track'):
			continue # we ignore these types of lines....
		else:
			parts = re.match(r"^(\w+)\s+(\d+)\s+(\d+)\s+(?:.+?)\s+(-?[\.\de-]+)", line)
			chrom = parts.group(1)
			start = int(parts.group(2))
			stop = int(parts.group(3))
			value = parts.group(4)
			output.write("%s\t%d\t%d\t%s\n" % (chrom, start, stop, value))
#end wig_to_bedgraph()

class Interval:
	def __init__(self, chrom, start, stop, score):
		self.chr = chrom
		self.start = start
		self.stop = stop
		self.score = score
	#end Interval.__init__()
	
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
			return Interval(self.chr, other.start, self.stop, 0)
	#end Interval.getOverlap()
	
	def __radd__(self, other):
		return other + self.score
	
	def __abs__(self):
		return abs(self.score)
#end class Interval

def convert_bdg_to_bw(bdg, chrsizes):
	
	with tempfile.NamedTemporaryFile() as temp_bdg:
		temp_bdg.write(bdg)
		temp_bdg.flush()
		with tempfile.NamedTemporaryFile() as temp_chrsizes:
			write_chrom_sizes(chrsizes, temp_chrsizes)
			temp_chrsizes.flush()
			with tempfile.NamedTemporaryFile() as temp_bw:
				p = subprocess.Popen(['bedGraphToBigWig', temp_bdg.name, temp_chrsizes.name, temp_bw.name], stderr=sys.stderr)
				p.communicate()
				temp_bw.seek(0)
				return temp_bw.read()
#end convert_bdg_to_bw()

def get_chrom_sizes(genome):
	if os.path.isfile(genome):
		return parse_chrom_sizes(genome)
	else:
		return fetch_chrom_sizes(genome)
#end get_chrom_sizes()

def write_chrom_sizes(chrom_sizes, outhandle):
	for chrom, size in chrom_sizes.iteritems():
		outhandle.write('%s\t%d\n' % (chrom, size))
#end write_chrom_sizes()

def parse_chrom_sizes(fname):
	sizes = {}
	with open(fname, 'r') as infile:
		reader = csv.reader(infile, delimiter='\t')
		for row in reader:
			sizes[row[0]] = int(row[1])
	return sizes
#end parse_chrom_sizes()

def fetch_chrom_sizes(genome_name):
	import subprocess
	infile = subprocess.check_output([	'mysql', 
										'--user=genome', 
										'--host=genome-mysql.cse.ucsc.edu', 
										'-ABN', '-D', genome_name, '-e', 'select chrom,size from chromInfo'
									], stderr=sys.stderr)

	sizes = {}
	rows = csv.reader(cStringIO.StringIO(infile), delimiter='\t')
	for row in rows:
		sizes[row[0]] = int(row[1])
	return sizes
#end fetch_chrom_sizes()

def parse_bedgraph(f, sort=True):
	intervals = []
	for line in f:
		match = re.search(r"(.+)\s+(\d+)\s+(\d+)\s+(-?[\d\.e-]+)", line)
		if match is not None:
			intervals.append(Interval(match.group(1), int(match.group(2)), int(match.group(3)), float(match.group(4))))
		else:
			sys.stderr.write('WARNING: Found non-conforming line: %s\n' % (line,))
	if sort:
		intervals.sort(key=lambda x:(x.chr, x.start, x.stop))
	return intervals
#end parse_bedgraph()

def write_bedgraph(intervals, outhandle):
	for interval in intervals:
		outhandle.write("%s\t%d\t%d\t%f\n" % (interval.chr, interval.start, interval.stop, interval.score))
#end write_bedgraph()

def repair_overlapping_segments(f, chroms, scorefunc):
	sys.stderr.write('Repairing overlapping regions\n')
	all_regions = parse_bedgraph(f)
	
	trees = {}
	for chrom in chroms:
		trees[chrom] = IntervalTree()
	
	for i in range(len(all_regions)):
		trees[all_regions[i].chr].addi(all_regions[i].start, all_regions[i].stop, all_regions[i].score)

	results = []
	for chrom in chroms:
		sys.stderr.write(' -> working on %s....\n' % (chrom,))
		trees[chrom].split_overlaps()
		trees[chrom].merge_equals(lambda l, v: l + [v], [])
	
		for i in trees[chrom]:
			results.append(Interval(chrom, i.begin, i.end, scorefunc(i.data)))
			
	results.sort(key=lambda x:(x.chr, x.start, x.stop))
	return results
#end repair_overlapping_segments()

def fill_complement(f, chroms, mode):
	"""fills the complement of the bedgraph with some value
	
	Args:
		f (string): b
	"""
	sys.stderr.write('Finding regions with missing data\n')
	
	if mode == 'zero':
		complement_value = 0
	else:
		complement_value = None
	
	all_regions = parse_bedgraph(f)
	genome_trees = {}
	for chrom in chroms:
		genome_trees[chrom] = IntervalTree()
		genome_trees[chrom].addi(0, chroms[chrom], complement_value)

	
	results = []
	for iv in all_regions:
		results.append(iv)
		genome_trees[iv.chr].chop(iv.start, iv.stop)
		
	for chrom in chroms:
		for iv in genome_trees[chrom]:
			results.append(Interval(chrom, iv.begin, iv.end, iv.data))	
			
	results.sort(key=lambda x:(x.chr, x.start, x.stop))
	return results
#end get_complement_fast()

