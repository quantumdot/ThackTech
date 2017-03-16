import os
import sys




class ProgressBar(object):


	def __init__(self, maxvalue, title="Progress", barlength=10, prettyDots=5, handle=None):
		self.barlength = barlength
		self.maxvalue = float(maxvalue)
		self.maxdots = prettyDots
		self.title = title
		self.laststatus = ""
		self.lastprogress = 0
		self.dotcount = 0
		self.maxstatuslength = 0
		self.isfirstprint = True
		if handle is None:
			self.handle = sys.stdout
		else:
			self.handle = handle

	def update(self, value, status=None, autoshow=True):
		self.lastprogress = value / self.maxvalue
		if status is not None:
			self.laststatus = status
		if autoshow:
			self.showProgress()

	def makeProgressBar(self):
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
			status += " Halt..."
		if progress >= 1:
			progress = 1
			posttext = "\r\n"
			status += " Done..."
		block = int(round(self.barlength*progress))
		dots = ("."*self.dotcount)+(" "*(self.maxdots-self.dotcount))
		text = self.title+"[{0}] {1:.2%} {2}".format( "#"*block + "-"*(self.barlength-block), progress, status+dots)
		#text += posttext
		return text

	def showProgress(self):
		if self.isfirstprint:
			pretext = ""
			self.isfirstprint = False
		else:
			pretext = "\033[2K\033[1G"
		self.handle.write(pretext+self.makeProgressBar())
		self.handle.flush()
#end class ProgressBar(object)