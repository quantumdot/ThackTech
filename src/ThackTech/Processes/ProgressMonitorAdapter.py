import os
import sys


class ProgressMonitorAdapter:
	def __init__(self, store_location):
		self.store_location = store_location
	def deleteAllProgress(self):
		pass
	def getStatus(self, task):
		pass
	def setStatus(self, task, status):
		pass
#end class ProgressMonitorAdapter()