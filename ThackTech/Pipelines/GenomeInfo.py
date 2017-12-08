import os
import glob
from ThackTech import filetools



class GenomeInfo(object):
	def __init__(self, name, gsize=0, chrsize=None):
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
			'BWAIndex': 	'*.fa.bwt',
			'Hisat2Index':  '*_tran.*.ht2',
			#'Hisat2Index':  'genome.*.ht2'
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
		if self.chrsize is None:
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
	@staticmethod
	def get_reference_genomes():
		if GenomeInfo.__known_references is None:
			GenomeInfo.__known_references = {}
			
			from ThackTech import conf
			genome_config = conf.get_config('genomes')
			for section in genome_config.sections():
				gi = GenomeInfo(section)
				
				#run the iGenomes discovery first.... possible to override with later directives.
				if genome_config.has_option(section, "goldenpath"):
					gi.try_discover(genome_config.get(section, "goldenpath"))
					
				options = genome_config.items(section)
				for oname, ovalue in options:
					if oname.startswith('index.'):
						idx_name = oname.split('.')[1]
						gi.add_index(idx_name, ovalue)
						
					elif oname.startswith('fasta.'):
						fa_name = oname.split('.')[1]
						if fa_name.lower() == 'genome':
							gi.wg_fasta = ovalue
						else:
							gi.chr_fasta[fa_name] = ovalue
					elif oname == 'size':
						gi.gsize = int(float(ovalue))
					else:
						setattr(gi, oname, ovalue)
				
				GenomeInfo.__known_references[gi.name] = gi
				
		return GenomeInfo.__known_references
	#end get_reference_genomes()
#end class GenomeInfo




