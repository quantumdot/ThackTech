import platform
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class StringTie(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='StringTie', short_description='Transcript assembly with StringTie.')
		super_args.update(**kwargs)
		super(StringTie, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('stringtie_path', str, 'stringtie', desc="Path to the StringTie executable."))
		
		self.add_parameter(ModuleParameter('min_abd', float, 0.1, desc="Sets the minimum isoform abundance of the predicted transcripts as a fraction of the most abundant transcript assembled at a given locus. See -f from stringtie."))
		self.add_parameter(ModuleParameter('min_len', int, 200, desc="Sets the minimum length allowed for the predicted transcripts. See -m from stringtie."))
		self.add_parameter(ModuleParameter('ballgown', bool, True, desc="Output Ballgown compatible files."))
		
		self.add_parameter(ModuleParameter('merge_mode', bool, False, desc="Run stringtie in transcript merge mode."))
		
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('alignments')
	#end __declare_resolvers()
	
	def tool_versions(self):
		return {
			'stringtie': self._call_output([self.get_parameter_value_as_string('stringtie_path'), '--version'], stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		cxt.log.write("Hello World from:\n{}\n".format(platform.uname()))
		cxt.log.flush()
		
		st_args = [
			self.get_parameter_value_as_string('stringtie_path'),
			'-p', str(self.processors),
			
			'-l', cxt.sample.name,
			'-f', self.get_parameter_value_as_string('min_abd'),
			'-m', self.get_parameter_value_as_string('min_len'),
		]
		
		if self.get_parameter_value('merge_mode'):
			st_args.append('--merge')
		
		else:
			
		
	#end run()
#end class StringTie