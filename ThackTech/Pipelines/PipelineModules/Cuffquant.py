import platform
from ThackTech.Pipelines import PipelineModule


class Cuffquant(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='Cuffquant', short_description='Quantify transcript expression with Cuffquant')
		super_args.update(**kwargs)
		super(Cuffquant, self).__init__(**super_args)
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