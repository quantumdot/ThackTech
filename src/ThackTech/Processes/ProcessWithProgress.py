import os
import sys
import subprocess
from ThackTech.Processes import ProgressStatusCodes



class ProcessWithProgress:
	
	def __init__(self, monitor, cmd, rollback_function=None, **kwargs):
		self.monitor = monitor
		self.process = None
		if hasattr(cmd, '__iter__'): #checks for all iterables except strings
			self.signature = " ".join(cmd)
		else:
			self.signature = cmd
		self.rollback_callback = rollback_function
		self.tryStartProcess(cmd, **kwargs)

	def tryStartProcess(self, cmd, **kwargs):
		if self.getStatus() == ProgressStatusCodes.INPROGRESS:
			self.doRollback()
		
		if self.getStatus() in [ ProgressStatusCodes.UNKNOWN, ProgressStatusCodes.NOTSTART ]:
			kwargs['preexec_fn'] = os.setpgrp
			if 'shell' in kwargs and kwargs['shell']:
				cmd = 'exec '+cmd
			self.process = subprocess.Popen(cmd, **kwargs)
			self.setStatus(ProgressStatusCodes.INPROGRESS)
				
	def doRollback(self):
		if self.rollback_callback is not None:
			self.rollback_callback()
		self.setStatus(ProgressStatusCodes.NOTSTART)
	
	def getStatus(self):
		return self.monitor.getStatus(self.signature)
	
	def setStatus(self, status):
		self.monitor.setStatus(self.signature, status)

	def isComplete(self):
		if self.getStatus() == ProgressStatusCodes.COMPLETED:
			return True
		return False
			
	def poll(self):
		if self.isComplete() is True:
			return 0
		ret = self.process.poll()
		if ret is not None:
			self.setStatus(ProgressStatusCodes.COMPLETED)
		return ret
		
	def communicate(self):
		if self.isComplete() is True:
			return 0
		(stdoutdata, stderrdata) = self.process.communicate()
		ret = self.process.returncode
		self.setStatus(ProgressStatusCodes.COMPLETED)
		return ret
#end class ProcessWithProgress			