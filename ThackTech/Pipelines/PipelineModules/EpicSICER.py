import os
import subprocess
import gzip
from ThackTech.Pipelines import PipelineModule, ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class EpicSICER(PipelineModule):
	
	def __init__(self, **kwargs):
		super(EpicSICER, self).__init__('EpicSICER', 'Peak Calling with Epic (SICER-based)', **kwargs)
		
		self.add_parameter(ModuleParameter('keep_duplicates', bool, False, desc="Keep reads mapping to the same position on the same strand within a library."))
		self.add_parameter(ModuleParameter('window_size', int, 200, desc="Size of the windows to scan the genome. WINDOW_SIZE is the smallest possible island."))
		self.add_parameter(ModuleParameter('gaps_allowed', int, 3, desc="Multiple of window size used to determine the gap size."))
		self.add_parameter(ModuleParameter('fragment_size', int, 3, desc="Size of the sequenced fragment. The center of the the fragment will be taken as half the fragment size."))
		self.add_parameter(ModuleParameter('fdr_cutoff', float, 1.0, desc="Remove all islands with an FDR below cutoff."))
		
		
		self._name_resolver('treatment')
		self._name_resolver('control')
	#end __init__()
	
	def tool_versions(self):
		return {
			'epic': self._call_output("epic --version 2>&1 | perl -ne 'if(m/epic ([\d\.-\w]+)/){ print $1; }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def load_modules(self, cxt):
		pass
	#end load_modules()
	
	def run(self, cxt):
		
		
		
		cmd = [
			'epic',
			'--number-cores', self.processors,
			'--genome', cxt.sample.genome.name,
			'--effective_genome_size', str(cxt.sample.genome.gsize),
			'--chromsizes', cxt.sample.genome.chrsize,
			
			'--window-size', self.get_parameter_value_as_string('window_size'),
			'--gaps-allowed', self.get_parameter_value_as_string('gaps_allowed'),
			
			'--false-discovery-rate-cutoff', self.get_parameter_value_as_string('fdr_cutoff'),
			
			
			'--store-matrix', _______MATRIX_FILENAME______,
			
		]
		if self.get_parameter_value('keep_duplicates'):
			cmd.append('--keep-duplicates')
		
		if cxt.sample.get_attribute('PE'):
			cmd.append('--paired-end')
		else:
			cmd.extend(['--fragment-size', self.get_parameter_value_as_string('fragment_size')])
		
		
		
		
	#end run()
#end class SamToBam