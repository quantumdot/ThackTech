import os
import sys
import time
from ThackTech import Common
from ThackTech.Processes import ProgressBar



class MultiStatusProgressItem:

	def __init__(self, name, status="", order=0):
		self.name = name
		self.order = order
		self.status = status
		self.progress = 0
		self.starttime = None
		self.endtime = None
		self.progressbar = ProgressBar(1, "", 10, 0)
		
	def start(self):
		self.starttime = time.time()
		return self

	def finish(self):
		self.endtime = time.time()
		return self
	
	def is_finished(self):
		return self.endtime is not None

	def elapsed_seconds(self):
		if self.starttime is None:
			return 0
		end = time.time() if self.endtime is None else self.endtime
		return end - self.starttime

	def elapsed_time(self):
		if self.starttime is None:
			return ""
		return Common.human_time_diff(self.starttime, time.time() if self.endtime is None else self.endtime)

	def update(self, progress=None, status=None):
		if progress is not None:
			self.progress = progress
		if status is not None:
			self.status = status
		self.progressbar.update(self.progress, self.elapsed_time(), False)
		return self

	def getData(self):
		self.update()
		return [self.name, self.status, self.progressbar.makeProgressBar()]
#end class MultiStatusProgressItem	