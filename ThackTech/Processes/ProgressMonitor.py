import os
import sys
from ThackTech.Processes import ProcessWithProgress


class ProgressMonitor:
	"""Provides an API for storing or quering the progress a script has made
	"""
	
	def __init__(self, store):
		self.store = store
	
	def deleteAllProgress(self):
		self.store.deleteAllProgress()

	def getStatus(self, task):
		return self.store.getStatus(task)
		
	def setStatus(self, task, status):
		self.store.setStatus(task, status)
		
	def PopenWithProgress(self, cmd, rollback, **kwargs):
		proc = ProcessWithProgress(self, cmd, rollback, **kwargs)
		return proc
	
	#def factory(self):
	#	return ProgressMonitor(self.store_location)
#end class ProgressMonitor