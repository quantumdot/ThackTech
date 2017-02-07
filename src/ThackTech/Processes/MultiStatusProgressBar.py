import os
import sys
import time
from tabulate import tabulate
from ThackTech import Common
from ThackTech.Processes import ProgressBar



class MultiStatusProgressBar(ProgressBar):
	

	def __init__(self, maxvalue, title=None, barlength=10, prettyDots=5, handle=None):
		super(MultiStatusProgressBar, self).__init__(maxvalue, title=title, barlength=barlength, prettyDots=prettyDots, handle=handle)
		self.lastmultistatus = {}
		self.lastheight = 0
		self.starttime = None
		self.endtime = None
		
	def update(self, value=None, globalstatus=None, taskstatuses=None):
		if taskstatuses is not None:
			self.lastmultistatus = taskstatuses
		super(MultiStatusProgressBar, self).update(value, globalstatus)
		return self
	
	def start(self):
		self.starttime = time.time()
		return self

	def finish(self):
		self.endtime = time.time()
		return self

	def elapsed_seconds(self):
		if self.starttime is None:
			return 0
		end = time.time() if self.endtime is None else self.endtime
		return end - self.starttime

	def elapsed_time(self):
		if self.starttime is None:
			return ""
		return Common.human_time_diff(self.starttime, time.time() if self.endtime is None else self.endtime)

	def showProgress(self):
		status = self.laststatus
		progress = self.lastprogress

		if self.dotcount >= self.maxdots:
			self.dotcount = 0
		else:
			self.dotcount += 1

		posttext = ""
		if isinstance(progress, int):
			progress = float(progress)
		if not isinstance(progress, float):
			progress = 0
			posttext = "\r\n"
			status = "error: progress var must be float!"
		if progress < 0:
			progress = 0
			posttext = "\r\n"
			status = "Halt..."
		if progress >= 1:
			progress = 1
			posttext = "\r\n"
			status = "Done..."

		block = int(round(self.barlength*progress))
		dots = ("."*self.dotcount)+(" "*(self.maxdots-self.dotcount))
			
		text = "\033[F\033[K" * self.lastheight
		text += self.formatTaskStatus()+"\n"
		text += self.title+": "
		text += "[{0}] {1:.2%} {2}\n".format('#'*block + "-"*(self.barlength-block), progress, status+dots)
		text += "Elapsed Time: %s\n" % self.elapsed_time()
		#text += posttext
		#text += self.formatTaskStatus()
		self.lastheight = (len(self.lastmultistatus)+4)
		self.handle.write(text)
		self.handle.flush()
	
	def formatTaskStatus(self):
		statuses = []
		if len(self.lastmultistatus) > 0:
			items = self.lastmultistatus.values()
			items.sort(key=lambda x:(x.order, x.name))
			for item in items:
				statuses.append(item.getData())
		return tabulate(statuses, headers=['Task', 'Status', 'Progress'], tablefmt="simple") 
#end class MultiStatusProgressBar(ProgressBar)			