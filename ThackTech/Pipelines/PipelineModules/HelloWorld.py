import platform
from ThackTech.Pipelines import PipelineModule


class HelloWorld(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='HelloWorld', short_description='Say "Hello World".')
		super_args.update(**kwargs)
		super(HelloWorld, self).__init__(**super_args)
	#end __init__()
	
	def run(self, cxt):
		cxt.log.write("Hello World from:\n%s\n" % (platform.uname(),))
		cxt.log.flush()
	#end run()
#end class HelloWorld