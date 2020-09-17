import os
import copy
import shutil
import filecmp
from ThackTech.Pipelines.Context import BaseModuleContext
from ThackTech import filetools

class FileContext(BaseModuleContext):
	"""Represents the context of a file within the pipeline system
    
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
	
	def __repr__(self):
		return "FileContext('%s', %d, '%s', '%s')" % (self.pipeline, self.step, self.module, self.role)
	#end __repr__()

#end class FileContext

class FileInfo(object):
	"""Represents a file within the pipeline framework
	
	Provides convenience methods for working with files, as well as providing support
	for tracing the generator of files. Also files can be associated with others
	through the use of companion files (i.e. a BAM file usually has a BAM index 
	file associated with it.
	
	FileInfo objects track the source of file through the use of three attributes:
	- generator: name of the module that produced the file
	- role: role of the file within the output (i.e. logfile, stats, etc.)
	- context: context that this file was generated under....... (NEED TO DEFINE THIS BETTER!!!)
	
	Important Nomenclature:
		- fullpath: /path/to/somefile.mult.ext
		- dirname:  /path/to/
		- basename: somefile.mult.ext
		- filename: somefile.mult
		- ext:      .ext
	
	"""
	def __init__(self, filepath, context=None, **attributes):
		"""Initializes a FileInfo object
		
		Parameters:
			filepath (str): Path of the file this object represents
		"""	
		#support for companion files
		from ThackTech.Pipelines import GLOBAL_MANAGER
		self.__companions = []
		
		#random extra data as needed.
		self.__context = context
		self.attributes = {}
		self.attributes.update(attributes)
		
		if isinstance(filepath, FileInfo):
			self.__fullpath = filepath.fullpath
			self.__companions.extend(filepath.companions)
		else:
			self.__fullpath = os.path.abspath(filepath)
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
		""" Shorthand for self.context
		"""
		return self.context
	
	@property
	def isfile(self):
		"""Tests if this is an existing regular file.
		
		see os.path.isfile() for more info
		"""
		return os.path.isfile(self.__fullpath)
	
	@property
	def exists(self):
		""" Tests if this file exists.
		
		see os.path.exists() for more info
		"""
		return os.path.exists(self.__fullpath)
	
	@property
	def fullpath(self):
		"""Gets the fully-qualified path to this file
		
		example:
		>>> f = FileInfo('/path/to/some/file.txt.gz')
		>>> f.fullpath
		>>> '/path/to/some/file.txt.gz'
		"""
		return os.path.abspath(self.__fullpath)
	
	@property
	def dirname(self):
		"""Gets the directory where this file is located
		Equivalent to os.path.dirname()
		
		example:
		>>> f = FileInfo('/path/to/some/file.txt.gz')
		>>> f.dirname
		>>> '/path/to/some'
		"""
		return os.path.dirname(self.__fullpath)
	
	@property
	def basename(self):
		"""Gets the name of this file, including the file extension
		Equivalent to os.path.basename()
		
		example:
		>>> f = FileInfo('/path/to/some/file.txt.gz')
		>>> f.basename
		>>> 'file.txt.gz'
		"""
		return os.path.basename(self.__fullpath)
		
	@property
	def filename(self):
		"""Gets the name of this file, without the terminal file extension
		Equivalent to os.path.splitext(os.path.basename())[0]
		
		example:
		>>> f = FileInfo('/path/to/some/file.txt.gz')
		>>> f.filename
		>>> 'file.txt'
		"""
		return os.path.splitext(self.basename)[0]
	
	@property
	def filename_strip_all_ext(self):
		"""Gets the name of this file, without ANY file extensions
		Equivalent to looping os.path.splitext(os.path.basename())[0]
		until no extensions are remaining
		
		example:
		>>> f = FileInfo('/path/to/some/file.txt.gz')
		>>> f.filename_strip_all_ext
		>>> 'file'
		"""
		return filetools.basename_noext(self.fullpath, complete=True)
		
	@property
	def ext(self):
		"""Gets the extension of this file
		Equivalent to os.path.splitext()[1]
		
		example:
		>>> f = FileInfo('/path/to/some/file.txt.gz')
		>>> f.ext
		>>> '.gz'
		"""
		return os.path.splitext(self.basename)[1]
	
	@property
	def companions(self):
		"""Gets a list of companion files associated with this file
		"""
		return self.__companions
	
	@property
	def has_companions(self):
		"""Returns True if this instance has companions, otherwise False
		"""
		return len(self.__companions) > 0
	
	def fullpath_with_ext(self, new_extension):
		"""Returns the fullpath of this file with the terminal file extension replaced by new_extension
		
		Parameters:
			new_extension: string file extension to replace terminal file extension with
		
		example:
		>>> f = FileInfo('/path/to/some/file.txt.gz')
		>>> f.fullpath_with_ext('bz2')
		>>> '/path/to/some/file.txt.bz2'
		"""
		return os.path.join(self.dirname, self.basename_with_ext(new_extension))
	
	def basename_with_ext(self, new_extension):
		"""Returns the basename of this file with the terminal file extension replaced by new_extension
		
		Parameters:
			new_extension: string file extension to replace terminal file extension with
		
		example:
		>>> f = FileInfo('/path/to/some/file.txt.gz')
		>>> f.basename_with_ext('bz2')
		>>> 'file.txt.bz2'
		"""
		return "{}.{}".format(self.filename, new_extension.lstrip('.'))
	
	def move(self, destfolder):
		""" Move the file that this instance represents, along with any companions
		
		Parameters:
			destfolder: string destination directory
			
		This method will move the file that this instance represents, including any
		companion files that are known. This instance will be updated to reflect the 
		new directory location of any file(s) moved. If a file with the same name 
		exists in the `destfolder`, that file is first removed in order to allow
		the move operation to proceed
		"""
		dest_file = os.path.join(destfolder, self.basename)
		
		#make sure the destination directory exists
		filetools.ensure_dir(destfolder)
		
		#check if the dest file already exists, if so remove it
		if os.path.exists(dest_file):
			os.remove(dest_file)
		
		#move the file
		shutil.move(self.fullpath, destfolder)
		self.__fullpath = dest_file
		
		#move any companions
		for companion in self.companions:
			companion.move(destfolder)
	
	def copy(self, destfolder):
		""" Copy the file that this instance represents, along with any companions
		
		Parameters:
			destfolder: string destination directory
			
		Returns:
			New FileInfo object representing the file copy.
			
		This method will copy the file that this instance represents, including any
		companion files that are known. A new FileInfo instance will be created and 
		returned to represent the new copy.
		"""
		shutil.copy2(self.fullpath, destfolder)
		fi = FileInfo(os.path.join(destfolder, self.basename), copy.deepcopy(self.context), **copy.deepcopy(self.attributes))
		for companion in self.companions:
			fi.companions.append(companion.copy(destfolder))
		return fi
		
	def __str__(self):
		return self.fullpath
	#end __str__()
	
	def __repr__(self):
		return "FileInfo('%s', %s)" % (self.fullpath, self.cxt.__repr__)
	#end __repr__()
	
	def __eq__(self, other):
		#lets see if other is (derived from) FileInfo 
		if not isinstance(self, other.__class__):
			return False
		
		#check if context is the same, if not consider it different
		if self.context != other.context:
			return False
		
		#check if both objs point to the same path
		if self.fullpath != other.fullpath:
			return False
		
		#check if attributes are the same
		if self.attributes != other.attributes:
			return False
		
		#check companion files
		if self.has_companions == other.has_companions: #both agree whether they have companions or not...
			
			if not self.has_companions: #both must not have companions
				return True
			else: #compare companions
				for c in self.companions:
					if c not in other.companions:
						return False #missing a companion
		else:
			return False #instance disagree on companions, so must be different	
		
		#with the above conditions satisfied, should only get here if everything is the same
		return True
	#end __eq__()

	def __ne__(self, other):
		return not self.__eq__(other)

#end class FileInfo







