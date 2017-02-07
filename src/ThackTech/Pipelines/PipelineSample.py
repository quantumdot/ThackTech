import os
import sys
from ThackTech.Pipelines import GenomeInfo, FileInfo, GLOBAL_MANAGER



class PipelineSample:

	def __init__(self, name, genome, dest, format):
		self.files = GLOBAL_MANAGER.dict()
		self.attr = GLOBAL_MANAGER.dict()
		self.name = name
		if isinstance(genome, GenomeInfo):
			self.genome = genome
		else:
			self.genome = GenomeInfo.get_reference_genomes()[genome]
		self.dest = os.path.abspath(dest)
		self.format = format.lower()
	#end __init__()
	
	def __getstate__(self):
		state = dict(self.__dict__)
		state['files'] = dict(self.files)
		state['attr'] = dict(self.attr)
		return state
	#end __getstate__()
	
	def __setstate__(self, state):
		self.files = GLOBAL_MANAGER.dict(state['files'])
		self.attr = GLOBAL_MANAGER.dict(state['attr'])
		self.__dict__.update(state)
	#end __setstate__()
	
	def add_file(self, group, label, path):
		if group not in self.files:
			self.files[group] = {}
		#necessary to swap the whole sub-dict when multiprocessing!
		tmp = self.files[group]
		tmp[label] = path
		self.files[group] = tmp
	#end add_file()
	
	def has_file_group(self, group):
		return group in self.files
	#end has_file_group()
	
	def has_file(self, group, label):
		if group in self.files:
			if label in self.files[group]:
				return True
		return False
	#end has_file()
	
	def get_file(self, group, label):
		if self.has_file(group, label):
			return self.files[group][label] 
		return None
	#end has_file()
	
	def get_file_group(self, group):
		if group in self.files:
			return self.files[group]
		return None
	#end get_file_group()
	
	def add_attribute(self, name, value):
		self.attr[name] = value
	#end add_attribute()
	
	def has_attribute(self, name):
		return name in self.attr
	#end has_attribute()
	
	def get_attribute(self, name):
		if self.has_attribute(name):
			return self.attr[name]
		return None
	#end get_attribute()
	
	def remove_attribute(self, name):
		if self.has_attribute(name):
			del self.attr[name]
		return None
	#end remove_attribute()
#end PipelineSample