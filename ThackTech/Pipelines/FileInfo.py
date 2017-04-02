import os
import shutil
import filecmp
from ThackTech.Pipelines import GLOBAL_MANAGER
from ThackTech.Pipelines.Context import BaseModuleContext

class FileContext(BaseModuleContext):
	"""Represents the context of a file within the pipelineing system
    
    """

	def __init__(self, pipeline_name, step_num, module_name, file_role):
		super(FileContext, self).__init__(pipeline_name, step_num, module_name)
		self.__role = file_role
	
	
	__origin_args = ("Pre-Pipeline", -1, "Origin")
	@staticmethod
	def from_origin(role):
		"""Creates a FileContext that represents the file was added during sample creation
		
		This is most often used to identify files that were input into the pipeline
		"""
		return FileContext(FileContext.__origin_args[0], FileContext.__origin_args[1], FileContext.__origin_args[2], role)
	#end from_origin
	
	@property
	def is_origin(self):
		return self.pipeline == FileContext.__origin_args[0] \
		   and self.step == FileContext.__origin_args[1] \
		   and self.module == FileContext.__origin_args[2]
	
	@staticmethod
	def from_module_context(context, role):
		"""Creates a FileContext from the supplied BaseModuleContext derivitive
		
		Parameters:
			context: (BaseModuleContext) Context to construct a new FileContext from
			role: (string) Role this file plays within the module context
			
		Returns:
			FileContext
		"""
		return FileContext(context.pipeline, context.step, context.module, role)
	#end from_module_context

	@property
	def role(self):
		"""The role of this file respect to the parent context (ModuleContext; ex: 'output_file_1')
		""" 
		return self.__role
	
	def __str__(self):
		return "%s_%s" % (super(FileContext, self).__str__(), self.role)

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
	def __init__(self, filepath, context=None, **attributes):
		"""Initializes a FileInfo object
		
		Args:
			filepath (str): Path of the file this object represents
		"""
		self.__fullpath = filepath
		self.__context = context
		
		#support for companion files
		self.__companions = GLOBAL_MANAGER.list()
		
		#random extra data as needed.
		self.attributes = {}
		self.attributes.update(attributes)
	#end __init__()
	
	def has_attribute(self, attribute):
		return attribute in self.attributes
	
	def has_attribute_value(self, attribute, value):
		if self.has_attribute(attribute):
			return self.attributes[attribute] == value
	
	def _set_path(self, filepath):
		"""Sets the internal file path of this FileInfo
		"""
		self.__fullpath = filepath
	
	@property
	def context(self):
		"""Gets the context of this file
		"""
		return self.__context
	@property
	def cxt(self):
		return self.context
	
	@property
	def fullpath(self):
		"""Gets the fully-qualified path to this file
		"""
		return os.path.abspath(self.__fullpath)
	
	@property
	def dirname(self):
		"""Gets the directory where this file is located
		Equivelent to os.path.dirname()
		"""
		return os.path.dirname(self.__fullpath)
	
	@property
	def basename(self):
		"""Gets the name of this file, including the file extension
		Equivelent to os.path.basename()
		"""
		return os.path.basename(self.__fullpath)
		
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
	
	@property
	def companions(self):
		"""Gets a list of companion files associated with this file
		"""
		return self.__companions
	
	def move(self, destfolder):
		shutil.move(self.fullpath, destfolder)
		self.__fullpath = os.path.join(destfolder, self.basename)
		for companion in self.companions:
			companion.move(destfolder)
	
	def copy(self, destfolder):
		shutil.copy2(self.fullpath, destfolder)
		for companion in self.companions:
			companion.move(destfolder)
		
	def __str__(self):
		return self.fullpath
	#end __str__()
	
	def __repr__(self):
		return "FileInfo(%s)" % (self.fullpath,)
	#end __repr__()
		
#end class FileInfo







