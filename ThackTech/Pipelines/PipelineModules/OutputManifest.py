import os
from ThackTech.Pipelines import PipelineModule



class OutputManifest(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='OutputManifest', short_description='Generate Output Manifest')
		super_args.update(**kwargs)
		super(OutputManifest, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		pass
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		pass
	#end __declare_resolvers()

	
	def run(self, cxt):
		cxt.log.write("-> Writing sample output manifest....\n")
		cxt.log.flush()
		cxt.sample.write_file_manifest(os.path.join(cxt.sample.dest, cxt.sample.name+'_output_manifest.tsv'))
	#end run()
#end class OutputManifest

