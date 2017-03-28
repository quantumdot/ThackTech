import platform
from ThackTech.Pipelines import PipelineModule


class HelloWorld(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'HelloWorld', 'Say "Hello World".')
	#end __init__()
	
	def run(self, cxt):
		cxt.log.write("Hello World from:\n%s\n" % (platform.uname(),))
		cxt.log.flush()
	#end run()
#end class HelloWorld