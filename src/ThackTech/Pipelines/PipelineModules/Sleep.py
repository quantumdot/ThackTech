import time
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class Sleep(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'Sleep', 'Take a nap.')
		
		self.add_parameter(ModuleParameter('sleep_time', int, 10, desc="Amount of time, in seconds, to sleep for."))
	#end __init__()

	
	def run(self, sample, logfile):
		logfile.write("About to sleep for %d seconds.\n" % (self.get_parameter_value('sleep_time'),))
		logfile.flush()
		time.sleep(self.get_parameter_value('sleep_time'))
	#end run()
#end class HelloWorld