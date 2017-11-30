import re
from ThackTech.Pipelines.PipelineModule import PipelineModule

class AnalysisPipeline(object):

	def __init__(self, name="Anonymous Pipeline"):
		self.name = name
		self.pipeline = []
		self.offset = None #or checkpoint name to start from
	#end __inti__()
	
# 	@property
# 	def num_steps(self):
# 		""" Gets the number of steps in this pipeline, including any offset
# 		"""
# 		return len(self.pipeline) + self.offset
	
	@property
	def safe_name(self):
		"""Gets a safe version of this pipeline's name
		
		The returned string should be safe to use in things like filenames etc.
		"""
		return re.sub('[^\w\-_\.]', '_', self.name)
	#end safe_name()
	
	def itersteps(self):
		i = 0
		started = self.offset == None
		steps_to_run = []
		for step in self.pipeline:
			if isinstance(step, AnalysisPipelineCheckpoint):
				if self.offset is not None and self.offset == step.name:
					started = True
				
			elif isinstance(step, PipelineModule):
				if started:
					steps_to_run.append(AnalysisPipelineStep(self.name, i, step))
				i += 1
		return steps_to_run
		#return [AnalysisPipelineStep(self.name, i, self.pipeline[i]) for i in range(len(self))]
	#end itersteps()
	
	def append_module(self, module, critical=None):
		"""Appends the module to the end of the pipeline
		"""
		if critical is not None:
			module.set_critical(bool(critical))
		self.pipeline.append(module)
	#end append_module()
	
	def insert_module(self, position, module, critical=None):
		"""Inserts the module at the specified index of the pipeline
		"""
		if critical is not None:
			module.set_critical(bool(critical))
		self.pipeline.insert(position, module)
	#end append_module()
	
	def documentation(self):
		"""Return a string that documents the pipeline
		"""
		hash_length = 40
		buff = "{}\nBegin Pipeline: {}\n{}\n".format('='*hash_length, self.name, '='*hash_length)
		
		for m in self.pipeline:
			buff += "\n{}\n".format(' | '*(hash_length/3))
			buff += "{}\n\n".format(' V '*(hash_length/3))
			buff += m.documentation()
		
		buff += "\n{}\n".format(' | '*(hash_length/3))
		buff += "{}\n\n".format(' V '*(hash_length/3))
		buff += "{}\nEnd Pipeline: {}\n{}\n".format('='*hash_length, self.name, '='*hash_length)
		return buff
	#end explain()
	
	def __len__(self):
		c = 0
		for s in self.pipeline:
			if isinstance(s, PipelineModule):
				c += 1
		return c
#end class AnalysisPipeline


class AnalysisPipelineStep(object):

	def __init__(self, pipeline_name, index, module_instance):
		self.__pipeline = pipeline_name
		self.__index = index
		self.__module = module_instance
	
	@property
	def pipeline(self):
		return self.__pipeline
	
	@property
	def index(self):
		return self.__index
	
	@property
	def module(self):
		return self.__module
	
	def __str__(self):
		return "{}_{}_{}".format(self.pipeline, self.step, self.module.name)
#end class AnalysisPipelineStep


class AnalysisPipelineCheckpoint(object):
	
	def __init__(self, checkpoint_name):
		self.name = checkpoint_name
		
	def documentation(self):
		"""Return a string that documents this checkpoint
		"""
		hash_length = 40
		buff = 'Pipeline Checkpoint\n'
		buff += '-' * hash_length
		buff += '\nCHECKPOINT NAME: {}\n'.format(self.name)
		buff += '-' * hash_length
		buff += '\n'
		return buff
#end class AnalysisPipelineCheckpoint




