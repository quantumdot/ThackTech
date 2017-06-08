import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext
from ThackTech import filetools


class EpicSICER(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='EpicSICER', short_description='Peak Calling with Epic (SICER-based)')
		super_args.update(**kwargs)
		super(EpicSICER, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('epic_path', str, 'epic', desc="path to EPIC"))
		self.add_parameter(ModuleParameter('keep_duplicates', bool, False, desc="Keep reads mapping to the same position on the same strand within a library."))
		self.add_parameter(ModuleParameter('window_size', int, 200, desc="Size of the windows to scan the genome. WINDOW_SIZE is the smallest possible island."))
		self.add_parameter(ModuleParameter('gaps_allowed', int, 3, desc="Multiple of window size used to determine the gap size."))
		self.add_parameter(ModuleParameter('fragment_size', int, 300, desc="Size of the sequenced fragment. The center of the the fragment will be taken as half the fragment size."))
		self.add_parameter(ModuleParameter('fdr_cutoff', float, 1.0, desc="Remove all islands with an FDR below cutoff."))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('treatments')
		self._name_resolver('controls')
	#end __declare_resolvers()
	
	def tool_versions(self):
		return {
			'epic': self._call_output(self.get_parameter_value('epic_path')+" --version 2>&1 | perl -ne 'if(m/epic ([\d\.-\w]+)/){ print $1; }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def load_modules(self, cxt):
		pass
	#end load_modules()
	
	def run(self, cxt):
		
		cmd = [
			self.get_parameter_value('epic_path'),
			'--number-cores', str(self.processors),
			'--genome', cxt.sample.genome.name,
			#'--effective_genome_size', str(cxt.sample.genome.gsize),
			#'--chromsizes', cxt.sample.genome.chrsize,
			
			'--window-size', self.get_parameter_value_as_string('window_size'),
			'--gaps-allowed', self.get_parameter_value_as_string('gaps_allowed'),
			'--false-discovery-rate-cutoff', self.get_parameter_value_as_string('fdr_cutoff'),
			'--sum-bigwig', cxt.sample.dest,
			'--bigwig', cxt.sample.dest,
			'--bed', os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.bed'),
			'--store-matrix', os.path.join(cxt.sample.dest, cxt.sample.name+'.matrix'),
		]
		if self.get_parameter_value('keep_duplicates'):
			cmd.append('--keep-duplicates')
		
		if cxt.sample.get_attribute('PE'):
			cmd.append('--paired-end')
		else:
			cmd.extend(['--fragment-size', self.get_parameter_value_as_string('fragment_size')])
		
		treatment_files = self.resolve_input('treatments', cxt)
		cmd.extend(['--treatment'] + [f.fullpath for f in treatment_files])
		
		control_files = self.resolve_input('controls', cxt)
		if len(control_files) > 0:
			cmd.extend(['--control'] + [f.fullpath for f in control_files])
		
		cxt.log.write("\t-> Performing peak calling with Epic......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(cmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		with open(os.path.join(cxt.sample.dest, cxt.sample.name+'_results.csv'), 'w') as results_out:
			self._run_subprocess(cmd, stderr=cxt.log, stdout=results_out)
		
		results = []
		results.append(FileInfo(os.path.join(cxt.sample.dest, cxt.sample.name+'_results.csv'), FileContext.from_module_context(cxt, 'results')))
		results.append(FileInfo(os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.bed'), FileContext.from_module_context(cxt, 'peaks')))
		results.append(FileInfo(os.path.join(cxt.sample.dest, cxt.sample.name+'.matrix'), FileContext.from_module_context(cxt, 'matrix')))
		results.append(FileInfo(os.path.join(cxt.sample.dest, 'chip_sum.bw'), FileContext.from_module_context(cxt, 'chip_sum_signal')))
		results.append(FileInfo(os.path.join(cxt.sample.dest, 'input_sum.bw'), FileContext.from_module_context(cxt, 'input_sum_signal')))
		
		for f in treatment_files:
			results.append(FileInfo(os.path.join(cxt.sample.dest, f.basename+'.bw'), FileContext.from_module_context(cxt, 'signal')))
		
		if len(control_files) > 0:
			for f in control_files:
				results.append(FileInfo(os.path.join(cxt.sample.dest, f.basename+'.bw'), FileContext.from_module_context(cxt, 'signal')))
		return results
	#end run()
#end class SamToBam