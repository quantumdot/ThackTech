import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class Cuffquant(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='Cuffquant', short_description='Quantify transcript expression with Cuffquant')
		super_args.update(**kwargs)
		super(Cuffquant, self).__init__(**super_args)
	#end __init__()
	
	def tool_versions(self):
		cq_path = self.get_parameter_value_as_string('cuffquant_path')
		return {
			'cuffquant': self._call_output(cq_path+" 2>&1 | perl -ne 'if(m/.*cuffquant.*v([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('cuffquant_path', str, 'cuffquant', desc="Path to the cuffquant executable."))
		
		ltype_opts = ['fr-unstranded', 'ff-firststrand', 'ff-secondstrand', 'ff-unstranded', 'fr-firststrand', 'fr-secondstrand', 'transfrags']
		self.add_parameter(ModuleParameter('library_type', str, ltype_opts[0], choices=ltype_opts, desc="Sets the library type. See --library-type from cuffquant."))
		
		
		
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		pass
	#end __declare_resolvers()
	
	def run(self, cxt):
		cq_args = [
			self.get_parameter_value_as_string('cuffquant_path'),
			'--quiet',
			'--num-threads', str(self.processors),
			'--output-dir', cxt.sample.dest,
			
		]
		
		
		
		
		
		
		
	#end run()
#end class HelloWorld