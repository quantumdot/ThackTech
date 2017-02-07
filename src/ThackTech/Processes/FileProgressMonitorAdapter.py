import os
import sys
import hashlib
import glob
from ThackTech.Processes import ProgressMonitor, ProgressMonitorAdapter, ProgressStatusCodes



class FileProgressMonitorAdapter(ProgressMonitorAdapter):
	
	def __init__(self, store_location, prefix=".", postfix=".progress"):
		ProgressMonitorAdapter.__init__(self, store_location)
		self.prefix = prefix
		self.postfix = postfix
	
	def __getfilepath(self, task):
		digest = hashlib.sha1(str(task).encode()).hexdigest()
		return os.path.join(self.store_location, self.prefix + digest + self.postfix)
	
	def deleteAllProgress(self):
		for f in glob.glob(os.path.join(self.store_location, self.prefix + '*' + self.postfix)):
			os.remove(f)

	def getStatus(self, task):
		try:
			with open(self.__getfilepath(task), 'r') as f:
				return int(f.read())
		except Exception:
			return ProgressStatusCodes.UNKNOWN
		
	def setStatus(self, task, status):
		with open(self.__getfilepath(task), 'w') as f:
			f.write(str(status))
#end class FileProgressMonitorAdapter()