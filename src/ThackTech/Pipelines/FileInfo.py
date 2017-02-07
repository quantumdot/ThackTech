import os
import sys
from ThackTech.Pipelines import GLOBAL_MANAGER

class FileInfo:
	"""Represents a file within the pipeline framework
	
	Provides convience methods for working with files, as well as providing support
	for tracing the generator of files. Also files can be associated with others
	through the use of companion files (i.e. a BAM file ususally has a BAM index 
	file associated with it.
	"""
	def __init__(self, filepath):
		"""Initializes a FileInfo object
		
		Args:
			filepath (str): Path of the file this object represents
		"""
		self._fullpath = filepath
		self._companions = GLOBAL_MANAGER.list()
		self._generator = ""
		
	
	@property
	def fullpath(self):
		return self._fullpath
	
	@property
	def dirname(self):
		return os.path.dirname(self._fullpath)
	
	@property
	def basename(self):
		return os.path.basename(self._fullpath)
		
	@property
	def filename(self):
		return os.path.splitext(self.basename)[0]
		
	@property
	def ext(self):
		return os.path.splitext(self.basename)[1]
	
	
	def has_companions(self):
		return len(self._companions) > 0
		
	def add_companion(self, role, filepath):
		self._companions.append(FileInfo(filepath))
#end class FileInfo