import os
import glob
from ThackTech import filetools



class GenomeInfo(object):
	def __init__(self, name, gsize, chrsize=None):
		"""Initalize a GenomeInfo object
		
		Parameters:
			name (str): genome name, typically UCSC naming convention
			gsize (int): Effective genome size
			chrsize (str): location of chromosome size info file
			
		"""
		self.name = name
		self.gsize = gsize
		self.chrsize = chrsize
		self.indexes = {}
		self.wg_fasta = None
		self.chr_fasta = {}
	#end __init__()
	
	def try_discover(self, basepath):
		'''given a base path, try to discover indexes and other reference genome data according to the illumina golden path layout'''
		#find indexes
		idx_types = {
			'BowtieIndex': 	'*.1.ebwt',
			'Bowtie2Index':	'*.1.bt2', 
			'BWAIndex': 	'*.fa.bwt'
		}
		for idx in idx_types.keys():
			idx_path = os.path.join(basepath, 'Sequence', idx)
			if os.path.exists(idx_path):
				matches = glob.glob(os.path.join(idx_path, idx_types[idx]))
				if len(matches) > 0:
					self.add_index(idx, os.path.join(idx_path, filetools.basename_noext(matches[0], True)))
		
		#chromosome fasta files
		for f in glob.glob(os.path.join(basepath, 'Sequence', 'Chromosomes', '*.fa')):
			self.chr_fasta[filetools.basename_noext(f)] = f
		
		#whole genome fasta
		wg_fasta_results = glob.glob(os.path.join(basepath, 'Sequence', 'WholeGenomeFasta', '*.fa'))
		if len(wg_fasta_results) > 0:
			self.wg_fasta = wg_fasta_results[0]
		
		#chromsizes
		if self.chrsize is not None:
			chrsize_results = glob.glob(os.path.join(basepath, 'Sequence', 'WholeGenomeFasta', 'chrom.sizes'))
			if len(chrsize_results) > 0:
				self.chrsize = chrsize_results[0]
	#end try_discover()
	
	def add_index(self, name, value):
		self.indexes[name] = value
	#end add_index
	
	def has_index(self, name):
		return name in self.indexes
	#end has_index()
	
	def get_index(self, name):
		if self.has_index(name):
			return self.indexes[name]
		return None
	#end get_index()
	
	__known_references = None
	__probe_paths = []
	@staticmethod
	def register_path(path):
		GenomeInfo.__probe_paths.append(path)
	#end register_path()
	
	@staticmethod
	def get_reference_genomes():
		if GenomeInfo.__known_references is None:
			GenomeInfo.__known_references = {}
			
# 			conanocal = []
# 			conanocal.append(GenomeInfo('hg18', '2.77e9'))
# 			conanocal.append(GenomeInfo('hg19', '2.79e9'))
# 			conanocal.append(GenomeInfo('mm9',  '1.91e9'))
# 			conanocal.append(GenomeInfo('mm10', '2.15e9'))
# 			
# 			for c in conanocal:
# 				for p in GenomeInfo.__probe_paths:
# 					#print "trying path %s for genome %s" % (os.path.join(p, c.name), c.name)
# 					if os.path.exists(os.path.join(p, c.name)):
# 						#print "Found!"
# 						c.try_discover(os.path.join(p, c.name))
# 						GenomeInfo.__known_references[c.name] = c
# 						break
			from ThackTech import conf
			genome_config = conf.get_config('genomes')
			for section in genome_config.sections():
				gi = GenomeInfo(section, int(genome_config.getfloat(section, "size")))
				if genome_config.has_option(section, "goldenpath"):
					gi.try_discover(genome_config.get(section, "goldenpath"))
				GenomeInfo.__known_references[gi.name] = gi
				
		return GenomeInfo.__known_references
	#end get_reference_genomes()
#end class GenomeInfo




# GenomeInfo.register_path('/mnt/ref/reference/Homo_sapiens/UCSC')
# GenomeInfo.register_path('/mnt/ref/reference/Mus_musculus/UCSC')
# GenomeInfo.register_path('/mnt/ref/reference/Rattus_norvegicus/UCSC')
# GenomeInfo.register_path('/mnt/ref/reference/Saccharomyces_cerevisiae/UCSC')
# GenomeInfo.register_path('/mnt/ref/reference/PhiX/UCSC')
# 
# GenomeInfo.register_path('/home/thackray/reference/Homo_sapiens/UCSC')
# GenomeInfo.register_path('/home/thackray/reference/Mus_musculus/UCSC')
# GenomeInfo.register_path('/home/thackray/reference/Rattus_norvegicus/UCSC')
# GenomeInfo.register_path('/home/thackray/reference/Saccharomyces_cerevisiae/UCSC')
# GenomeInfo.register_path('/home/thackray/reference/PhiX/UCSC')




