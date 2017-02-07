import os
import platform
import sys
import time

import ThackTech.Common as Common
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class HelloWorld(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'HelloWorld', 'Say "Hello World" and take a nap.')
		
		self.add_parameter(ModuleParameter('sleep_time', int, 10, desc="Amount of time, in seconds, to sleep after saying hello."))
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, sample, logfile):
		logfile.write("Hello World from:\n%s\n" % (platform.uname(),))
		logfile.flush()
		time.sleep(self.get_parameter_value('sleep_time'))
	#end run()
#end class HelloWorld