import platform
from ThackTech.Pipelines import PipelineModule


class HelloWorld(PipelineModule):
	
	def __init__(self, **kwargs):
		super(HelloWorld, self).__init__('HelloWorld', 'Say "Hello World".', **kwargs)
	#end __init__()
	
	def run(self, cxt):
		cxt.log.write("Hello World from:\n%s\n" % (platform.uname(),))
		cxt.log.flush()
	#end run()
#end class HelloWorld