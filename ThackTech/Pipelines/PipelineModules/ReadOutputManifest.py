import ast
import os
from collections import OrderedDict
from ThackTech.Pipelines import PipelineModule
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class ReadOutputManifest(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='ReadOutputManifest', short_description='Read Output Manifest')
		super_args.update(**kwargs)
		super(ReadOutputManifest, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		pass
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		pass
	#end __declare_resolvers()
	
	def run(self, cxt):
		cxt.log.write("-> Reading output manifest for sample {}....\n".format(cxt.sample.name))
		cxt.log.flush()
		cxt.sample.read_file_manifest(os.path.join(cxt.sample.dest, cxt.sample.name+'_output_manifest.tsv'))
	#end run()
#end class ReadOutputManifest