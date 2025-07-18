import time
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class Sleep(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='Sleep', short_description='Take a nap.')
		super_args.update(**kwargs)
		super(Sleep, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('sleep_time', int, 10, desc="Amount of time, in seconds, to sleep for."))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		pass
	#end __declare_resolvers()

	
	def run(self, cxt):
		cxt.log.write("About to sleep for %d seconds.\n" % (self.get_parameter_value('sleep_time'),))
		cxt.log.flush()
		time.sleep(self.get_parameter_value('sleep_time'))
	#end run()
#end class HelloWorld