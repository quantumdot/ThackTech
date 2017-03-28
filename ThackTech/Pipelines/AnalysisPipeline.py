import re
from ThackTech.Pipelines.Context import ModuleContext

class AnalysisPipeline(object):

	def __init__(self, name="Anonymous Pipeline"):
		self.name = name
		self.pipeline = []
	#end __inti__()
	
	@property
	def safe_name(self):
		"""Gets a safe version of this pipeline's name
		
		The returned string should be safe to use in things like filenames etc.
		"""
		return re.sub('[^\w\-_\.]', '_', self.name)
	#end safe_name()
	
	def itersteps(self):
		return [ModuleContext(self.name, i, self.pipeline[i]) for i in range(len(self))]
	#end itersteps()
	
	def append_module(self, module, critical=False):
		"""Appends the module to the end of the pipeline
		"""
		module.set_critical(critical)
		self.pipeline.append(module)
	#end append_module()
	
	def insert_module(self, position, module, critical=False):
		"""Inserts the module at the specified index of the pipeline
		"""
		module.set_critical(critical)
		self.pipeline.insert(position, module)
	#end append_module()
	
	def documentation(self):
		"""Return a string that documents the pipeline
		"""
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
	
	def __len__(self):
		return len(self.pipeline)
#end class AnalysisPipeline
