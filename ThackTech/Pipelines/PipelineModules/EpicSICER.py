import os
import subprocess
import gzip
from ThackTech.Pipelines import PipelineModule, ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext
from ThackTech import filetools
from IPython.utils.io import stdout


class EpicSICER(PipelineModule):
	
	def __init__(self, **kwargs):
		super(EpicSICER, self).__init__('EpicSICER', 'Peak Calling with Epic (SICER-based)', **kwargs)
		
		self.add_parameter(ModuleParameter('keep_duplicates', bool, False, desc="Keep reads mapping to the same position on the same strand within a library."))
		self.add_parameter(ModuleParameter('window_size', int, 200, desc="Size of the windows to scan the genome. WINDOW_SIZE is the smallest possible island."))
		self.add_parameter(ModuleParameter('gaps_allowed', int, 3, desc="Multiple of window size used to determine the gap size."))
		self.add_parameter(ModuleParameter('fragment_size', int, 300, desc="Size of the sequenced fragment. The center of the the fragment will be taken as half the fragment size."))
		self.add_parameter(ModuleParameter('fdr_cutoff', float, 1.0, desc="Remove all islands with an FDR below cutoff."))
		
		
		self._name_resolver('treatments')
		self._name_resolver('controls')
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
		
		filetools.ensure_dir(os.path.join(cxt.sample.dest, 'signal'))
		filetools.ensure_dir(os.path.join(cxt.sample.dest, 'sum_signal'))
		
		
		cmd = [
			'epic',
			'--number-cores', self.processors,
			'--genome', cxt.sample.genome.name,
			#'--effective_genome_size', str(cxt.sample.genome.gsize),
			#'--chromsizes', cxt.sample.genome.chrsize,
			
			'--window-size', self.get_parameter_value_as_string('window_size'),
			'--gaps-allowed', self.get_parameter_value_as_string('gaps_allowed'),
			'--false-discovery-rate-cutoff', self.get_parameter_value_as_string('fdr_cutoff'),
			'--sum-bigwig', os.path.join(cxt.sample.dest, 'sum_signal'),
			'--bigwig', os.path.join(cxt.sample.dest, 'signal'),
			'--bed', os.path.join(cxt.sample.dest, 'bed'),
			'--store-matrix', os.path.join(cxt.sample.dest, 'matrix', cxt.sample.name+'.matrix'),
		]
		if self.get_parameter_value('keep_duplicates'):
			cmd.append('--keep-duplicates')
		
		if cxt.sample.get_attribute('PE'):
			cmd.append('--paired-end')
		else:
			cmd.extend(['--fragment-size', self.get_parameter_value_as_string('fragment_size')])
		
		treatment_files = self.resolve_input('treatments', cxt)
		cmd.extend(['--treatment']+treatment_files)
		
		control_files = self.resolve_input('controls', cxt)
		if len(control_files) > 0:
			cmd.extend(['--control']+control_files)
		
		with open(os.path.join(cxt.sample.dest, cxt.sample.name+'_results.csv')) as results_out:
			self._run_subprocess(cmd, stderr=cxt.log, stdout=results_out)
			
		
	#end run()
#end class SamToBam