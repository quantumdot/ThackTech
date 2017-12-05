import os
import platform
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class GffCompare(PipelineModule):
	'''
	STATUS: EXPERIMENTAL - Not Fully Implemented!
	'''
	
	def __init__(self, **kwargs):
		super_args = dict(name='GffCompare', short_description='Compare GFF files.')
		super_args.update(**kwargs)
		super(GffCompare, self).__init__(**super_args)
	#end __init__()
	
	def tool_versions(self):
		p = self.get_parameter_value_as_string('gffcompare_path')
		return {
			'gffcompare': self._call_output(p+" --version 2>&1 | perl -ne 'if(m/.*gffcompare.*v([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('gffcompare_path', str, 'gffcompare', desc="Path to the gffcompare executable."))
		
	#end __declare_parameters()
	
	def run(self, cxt):		
		pass
		
	#end run()

#end class GffCompare
