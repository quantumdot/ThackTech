import re

class AnalysisPipeline:

	def __init__(self, name="Anonymous Pipeline"):
		self.name = name
		self.pipeline = []
		self.tasks_statuses = None
	#end __inti__()
	
	def safe_name(self):
		return re.sub('[^\w\-_\.]', '_', self.name)
	#end safe_name()
	
	def append_module(self, module, critical=False):
		module.set_critical(critical)
		self.pipeline.append(module)
	#end append_module()
	
	def insert_module(self, position, module, critical=False):
		module.set_critical(critical)
		self.pipeline.insert(position, module)
	#end append_module()
	
	def documentation(self):
		hash_length = 40
		buff = "%s\nBegin Pipeline: %s\n%s\n" % ('='*hash_length, self.name, '='*hash_length)
		
		for m in self.pipeline:
			buff += "\n%s\n" % (u' | '*(hash_length/3),)
			buff += "%s\n\n" % (u' V '*(hash_length/3),)
			buff += m.documentation()
		
		buff += "\n%s\n" % (u' | '*(hash_length/3),)
		buff += "%s\n\n" % (u' V '*(hash_length/3),)
		buff += "%s\nEnd Pipeline: %s\n%s\n" % ('='*hash_length, self.name, '='*hash_length)
		return buff
	#end explain()
#end class AnalysisPipeline
