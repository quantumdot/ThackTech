import os
import subprocess
import gzip
from ThackTech.Pipelines import PipelineModule, ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class Conditional(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='Conditional', short_description='Run modules depending on some condition')
		super_args.update(**kwargs)
		super(Conditional, self).__init__(**super_args)
		
		self.true_children = []
		self.false_children = []
		
		self._name_resolver('condition')
	#end __init__()
	
	def tool_versions(self):
		return {}
	#end tool_versions()
	
	def load_modules(self, cxt):
		pass
	#end load_modules()
	
	def run(self, cxt):
		cxt.log.write("\t-> Converting BAM to BED...\n")
		
		
		
	#end run()
#end class SamToBam