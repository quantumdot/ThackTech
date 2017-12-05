import platform
from ThackTech.Pipelines import PipelineModule


class Cuffmerge(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='Cuffmerge', short_description='Merge Transcript assemblies with Cuffmerge')
		super_args.update(**kwargs)
		super(Cuffmerge, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		pass
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		pass
	#end __declare_resolvers()
	
	def run(self, cxt):
		cxt.log.write("Hello World from:\n{}\n".format(platform.uname()))
		cxt.log.flush()
	#end run()
#end class HelloWorld