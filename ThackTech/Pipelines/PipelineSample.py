import os
import ast
from collections import OrderedDict
from ThackTech.Pipelines import GenomeInfo, FileInfo, FileContext, GLOBAL_MANAGER



class PipelineSample(object):
	"""A PipelineSample represents an individual sample that moves through a pipeline.
	
	A PipelineSample serves as a container for sample data (i.e. state) as it traverses 
	a pipeline. samples contain data that identifies the sample, files associated with this
	sample, data regarding the reference genome this sample uses, and a dict object for storing
	arbitrary information.
	
	The file and attribute tables are managed by a multiprocessing manager to enable
	multiprocessing support. These should not be accessed directly in favor of the method
	wrappers.
	
	"""
	def __init__(self, name, genome, dest):
		"""Initializes a PipelineSample object
		
		Parameters:
			name:	(string)				Name of this sample
			genome:	(GenomeInfo | string)	String representing known reference genome or directly a GenomeInfo object
			dest:	(string)				Directory path that serves as the destination for results
		"""
		self.files = GLOBAL_MANAGER.list()
		self.attr = GLOBAL_MANAGER.dict()
		self.name = name
		if isinstance(genome, GenomeInfo):
			self.genome = genome
		else:
			self.genome = GenomeInfo.get_reference_genomes()[genome]
		self.dest = os.path.abspath(dest)
	#end __init__()
	
	def __getstate__(self):
		"""implemented for pickling support
		"""
		state = dict(self.__dict__)
		state['files'] = list(self.files)
		state['attr'] = dict(self.attr)
		return state
	#end __getstate__()
	
	def __setstate__(self, state):
		"""implemented for pickling support
		"""
		self.files = GLOBAL_MANAGER.list(state['files'])
		self.attr = GLOBAL_MANAGER.dict(state['attr'])
		self.__dict__.update(state)
	#end __setstate__()
	
	
	def set_attribute(self, name, value):
		"""Set an attribute and its value, overwriting any previous value if set
		
		Parameters:
			name:	(string) Name of the attribute
			value:	(??????) Value of the attribute
		"""
		self.attr[name] = value
	#end add_attribute()
	
	def has_attribute(self, name):
		"""Return a boolean indicating if name is in the attribute dictionary
		
		Parameters:
			name:	(string)	Name of the attribute to check
		
		Returns:
			(bool) True if name exists as a key in the attribute dictionary
		"""
		return name in self.attr
	#end has_attribute()
	
	def get_attribute(self, name):
		"""Gets the value of the attribute with name name, or None if name does not exist
		
		Parameters:
			name:	(string)	Name of the attribute to get
		
		Returns:
			(None | value) None if attribute does not exist, otherwise arbitrary value
		"""
		if self.has_attribute(name):
			return self.attr[name]
		return None
	#end get_attribute()
	
	def remove_attribute(self, name):
		"""deletes the attribute with name name
		
		Parameters:
			name:	(string)	Name of the attribute to delete
		"""
		if self.has_attribute(name):
			del self.attr[name]
	#end remove_attribute()
	
	
	
	
	def add_file(self, fileinfo, allow_duplicates=False):
		"""Adds a file to this sample
		
		Parameters:
			fileinfo:	String path to the file, or (preferably) a FileInfo object
		"""
		if not isinstance(fileinfo, FileInfo):
			fileinfo = FileInfo(fileinfo)
			
		if not allow_duplicates and fileinfo in self.files:
			return #duplicates not allowed
		
		self.files.append(fileinfo)
	#end add_file()
	
	def clear_files(self):
		while len(self.files) > 0:
			self.files.pop()
	
	def remove_file(self, fileinfo):
		self.files.remove(fileinfo)
	
	def find_files(self, predicate):
		"""Find files that satisfy a predicate
		
		Find and return a list of FileInfo objects 
		that satisfy the supplied predicate. The FileInfo
		object is supplied to the predicate as the sole argument.
		>>> sample.find_files(lambda f: f.ext == 'bam')
		
		Parameters:
			predicate: Callable that takes a FileInfo as its sole parameter
			
		Returns:
			List of FileInfo objects that satisfy the supplied predicate
		"""
		return list(filter(predicate, self.files))
	#end find_files()
	
	def write_file_manifest(self, path):
		
		def write_manifest_line(index, parent, f, out):
			tpl = "{id}\t{parent}\t{pipeline}\t{step}\t{module}\t{role}\t{path}\t{attributes}\n"
			out.write(tpl.format(id=index, parent=parent, pipeline=f.cxt.pipeline, step=f.cxt.step,
								 module=f.cxt.module, role=f.cxt.role, path=f.fullpath,
								 attributes=';'.join("%s=%r" % (key,val) for (key,val) in f.attributes.iteritems())))
		#end inner write_manifest_line()
		
		with open(path, 'w') as output_manifest:
			output_manifest.write('id\tparent\tpipeline\tstep\tmodule\trole\tpath\tattributes\n')
			
			i = 0
			for f in self.files:
				write_manifest_line(i, "None", f, output_manifest)
				
				if len(f.companions) > 0:
					j=1
					for fc in f.companions:
						write_manifest_line(i+j, i, fc, output_manifest)
						j += 1
					i += j
				
				i += 1
	#end write_file_manifest()
	
	def read_file_manifest(self, path):
		count = 0
		man_files = OrderedDict()
		with open(path, 'r') as output_manifest:
			for line in output_manifest:
				line = line.strip()
				if line == '':
					continue
				if line.startswith('id'):
					continue #header line
				parts = line.split('\t')
				if len(parts) > 0:
					c = FileContext(parts[2], int(parts[3]), parts[4], parts[5])
					f = FileInfo(parts[6], c)
					if len(parts) > 7 and parts[7].strip() != '':
						tuples = [item.split("=") for item in parts[7].split(";")]
						for t in tuples:
							f.attributes[t[0]] = ast.literal_eval(t[1])
					
					if parts[1] == 'None':
						#this is a primary file
						man_files[parts[0]] = f
					else:
						#this is a companion file
						man_files[parts[1]].companions.append(f)

					count += 1
		
		self.clear_files() #clear so we don't get duplicates
		for f in man_files.values():
			self.add_file(f)
	#end read_file_manifest()
	
	#===========================================================================
	# DEPRECATED
	#===========================================================================
	# def has_file_group(self, group):
	# 	return group in self.files
	# #end has_file_group()
	# 
	# def has_file(self, group, label):
	# 	if group in self.files:
	# 		if label in self.files[group]:
	# 			return True
	# 	return False
	# #end has_file()
	# 
	# def get_file(self, group, label):
	# 	if self.has_file(group, label):
	# 		return self.files[group][label] 
	# 	return None
	# #end has_file()
	# 
	# def get_file_group(self, group):
	# 	if group in self.files:
	# 		return self.files[group]
	# 	return None
	# #end get_file_group()
	#===========================================================================
#end PipelineSample