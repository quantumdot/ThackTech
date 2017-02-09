import os
import sys
from ThackTech.Pipelines import GLOBAL_MANAGER

class ModuleContext(object):

	def __init__(self, pipeline_name, step_num, module_name):
		self.__pipeline = pipeline_name
		self.__step = step_num
		self.__module = module_name
	
	@property
	def pipeline(self):
		return self.__pipeline
	
	@property
	def step(self):
		return self.__step
	
	@property
	def module(self):
		return self.__module
	
#end class ModuleContext

class ModuleRunContext(ModuleContext):

	def __init__(self, pipeline_name, step_num, module_name, logfile, sample):
		super(FileContext, self).__init__(pipeline_name, step_num, module_name)
		self.__logfile = logfile
		self.__sample = sample
		
	@property
	def sample(self):
		return self.__sample
	
	@property
	def logfile(self):
		return self.__logfile
	
#end class ModuleRunContext

class FileContext(ModuleContext):
	"""Represents the context of a file within the pipelineing system
	
	"""
	def __init__(self, pipeline_name, step_num, module_name, file_role):
		super(FileContext, self).__init__(pipeline_name, step_num, module_name)
		self.__role = file_role
		
	@property
	def role(self):
		return self.__role
	
#end class FileContext

class FileInfo(object):
	"""Represents a file within the pipeline framework
	
	Provides convience methods for working with files, as well as providing support
	for tracing the generator of files. Also files can be associated with others
	through the use of companion files (i.e. a BAM file ususally has a BAM index 
	file associated with it.
	
	FileInfo objects track the source of file through the use of three arrtibutes:
	- generator: name of the module that produced the file
	- role: role of the file within the output (i.e. logfile, stats, etc.)
	- context: context that this file was generated under....... (NEED TO DEFINE THIS BETTER!!!)
	
	"""
	def __init__(self, filepath):
		"""Initializes a FileInfo object
		
		Args:
			filepath (str): Path of the file this object represents
		"""
		self._fullpath = filepath
		
		#contextual information
		self._generator = ""
		self._role = ""
		self.context = ""
		
		#support for companion files
		self._parent = None
		self._companions = GLOBAL_MANAGER.list()
	#end __init__()
	
	@property
	def fullpath(self):
		"""Gets the fully-qualified path to this file
		"""
		return self._fullpath
	
	@property
	def dirname(self):
		"""Gets the directory where this file is located
		Equivelent to os.path.dirname()
		"""
		return os.path.dirname(self._fullpath)
	
	@property
	def basename(self):
		"""Gets the name of this file, including the file extension
		Equivelent to os.path.basename()
		"""
		return os.path.basename(self._fullpath)
		
	@property
	def filename(self):
		"""Gets the name of this file, without the file extension
		Equivelent to os.path.splitext(os.path.basename())[0]
		"""
		return os.path.splitext(self.basename)[0]
		
	@property
	def ext(self):
		"""Gets the extension of this file
		Equivelent to os.path.splitext()[1]
		"""
		return os.path.splitext(self.basename)[1]
	
	
	def has_companions(self):
		return len(self._companions) > 0
		
	def add_companion(self, role, filepath):
		self._companions.append(FileInfo(filepath))
		
		
		
		
#end class FileInfo







