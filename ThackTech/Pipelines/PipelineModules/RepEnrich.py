import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter, GenomeInfo


def __register_repenrich_index():
	genomes = GenomeInfo.get_reference_genomes()
	genomes['hg19'].add_index('repenrich_setup', '/mnt/ref/RepEnrich_data/setup_hg19')
	genomes['hg19'].add_index('repenrich_annot', '/mnt/ref/repeatmasker/hg19_repeatmasker.4.0.5.txt.fix')

	genomes['mm9'].add_index('repenrich_setup', '/mnt/ref/RepEnrich_data/setup_mm9')
	genomes['mm9'].add_index('repenrich_annot', '/mnt/ref/repeatmasker/mm9_repeatmasker.txt')
#end __register_repenrich_index()
__register_repenrich_index()



class RepEnrich(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='RepEnrich', short_description='Estimate Repetitive Element Enrichment')
		super_args.update(**kwargs)
		super(RepEnrich, self).__init__(**super_args)
	#end __init__()
	
	def __declare_parameters(self):
		self.add_parameter(ModuleParameter('repenrich_path', str, '/home/josh/scripts/RepEnrich3.py', desc="Location of RepEnrich script."))
		self.add_parameter(ModuleParameter('cleanup', bool, False, desc="Indicates if cleanup should be performed on intermediary files."))
		#self.add_parameter(ModuleParameter('cleanup', bool, False, desc="Indicates if cleanup should be performed on intermediary files."))
	#end __declare_parameters()
	
	def __declare_resolvers(self):
		self._name_resolver('multimap_fastq')
		self._name_resolver('unique_bam')
		self.set_resolver('annotation_file', lambda cxt: cxt.sample.genome.get_index('repenrich_annot'))
		self.set_resolver('setup_folder', lambda cxt: cxt.sample.genome.get_index('repenrich_setup'))
	#end __declare_resolvers()

	def supported_types(self):
		return None
	#end supported_types()
	
	def show_version(self, handle=None, fancy=True):
		pass
		# handle.write('\nBowtie version is:\n------------------------------------------\n')
		# handle.flush()
		# proc = subprocess.Popen('bowtie --version', shell=True)
		# proc.communicate()
		# handle.write('------------------------------------------\n\n')
		# handle.flush()
	#end show_version()
	
	def run(self, cxt):
		cxt.log.write("\t-> Preparing RepEnrich....\n")
		cxt.log.flush()
		
		multimap_fastq = self.resolve_input('multimap_fastq', cxt.sample)
		if isinstance(multimap_fastq, str):
			multimap_fastq = [multimap_fastq]
			
		repenrich_cmd = [
			'python', 		self.get_parameter_value_as_string('repenrich_path'),
			'--noshowprogress',
			'--progressdir', cxt.sample.dest,
			'--workers', 	str(self.processors),
			'--annotation',	self.resolve_input('annotation_file', cxt.sample),
			'--dest',		cxt.sample.dest,
			'--name',		cxt.sample.name,
			'--setup',		self.resolve_input('setup_folder', cxt.sample),
			'--bam',		self.resolve_input('unique_bam', cxt.sample)
		] + ['--fastq'] + multimap_fastq

		#OK, we now have all the arguments setup, lets actually run bowtie
		cxt.log.write("\t-> RepEnrich: Estimate Repetitive Element Enrichment......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(repenrich_cmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		proc = subprocess.Popen(repenrich_cmd, stderr=subprocess.STDOUT, stdout=cxt.log)
		proc.communicate()
		
		return {
			
		}
	#end run()
#end class BowtieAlign