import sqlite3
import subprocess
import sys
import os
from Queue import Queue, Empty
from threading import Thread, Event
import multiprocessing
import time
import Common
import hashlib
import glob

class ProgressMonitor:
	'Provides an API for storing or quering the progress a script has made'
	
	class StatusCodes:
		UNKNOWN = 0
		NOTSTART = 1
		INPROGRESS = 2
		COMPLETED = 3
	
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

class ProgressMonitorAdapter:
	def __init__(self, store_location):
		self.store_location = store_location
	def deleteAllProgress(self):
		pass
	def getStatus(self, task):
		pass
	def setStatus(self, task, status):
		pass
#end class FileProgressMonitorAdapter()


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
			return ProgressMonitor.StatusCodes.UNKNOWN
		
	def setStatus(self, task, status):
		with open(self.__getfilepath(task), 'w') as f:
			f.write(str(status))
#end class FileProgressMonitorAdapter()
		
		
class SQLiteProgressMonitorAdapter(ProgressMonitorAdapter):
	
	def __init__(self, store_location, table_name='progress', task_col='task', stat_col='status'):
		ProgressMonitorAdapter.__init__(self, store_location)
		self.table_name = table_name
		self.task_col = task_col
		self.stat_col = stat_col
		self.store = sqlite3.connect(self.store_location, timeout=300)
		self._initTable()
		
	def _initTable(self):
		sql = 'CREATE TABLE IF NOT EXISTS '+self.table_name+' ('+self.stat_col+' INTEGER, '+self.task_col+' TEXT PRIMARY KEY)'
		c = self.store.cursor()
		c.execute(sql)
		self.store.commit()
	
	def deleteAllProgress(self):
		sql = 'DROP TABLE '+self.table_name
		c = self.store.cursor()
		c.execute(sql)
		self.store.commit()
		self.initTable()

	def getStatus(self, task):
		c = self.store.cursor()
		c.execute('SELECT '+self.stat_col+' FROM '+self.table_name+' WHERE '+self.task_col+'=?', (task,))
		r = c.fetchone()
		if r is None:
			return ProgressMonitor.StatusCodes.UNKNOWN
		return r[0]
		
	def setStatus(self, task, status):
		c = self.store.cursor()
		c.execute('INSERT OR REPLACE INTO '+self.table_name+' ('+self.task_col+', '+self.stat_col+') VALUES (?, ?)', (task, status))
		self.store.commit()
#end class SQLiteProgressMonitorAdapter





class ProcessTicket:
	
	def __init__(self, signature, status=ProgressMonitor.StatusCodes.UNKNOWN):
		self.times = {}
		self.signature = signature
		self.status = status
#end class ProcessTicket

	
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
		if self.getStatus() == ProgressMonitor.StatusCodes.INPROGRESS:
			self.doRollback()
		
		if self.getStatus() in [ ProgressMonitor.StatusCodes.UNKNOWN, ProgressMonitor.StatusCodes.NOTSTART ]:
			kwargs['preexec_fn'] = os.setpgrp
			if 'shell' in kwargs and kwargs['shell']:
				cmd = 'exec '+cmd
			self.process = subprocess.Popen(cmd, **kwargs)
			self.setStatus(ProgressMonitor.StatusCodes.INPROGRESS)
				
	def doRollback(self):
		if self.rollback_callback is not None:
			self.rollback_callback()
		self.setStatus(ProgressMonitor.StatusCodes.NOTSTART)
	
	def getStatus(self):
		return self.monitor.getStatus(self.signature)
	
	def setStatus(self, status):
		self.monitor.setStatus(self.signature, status)

	def isComplete(self):
		if self.getStatus() == ProgressMonitor.StatusCodes.COMPLETED:
			return True
		return False
			
	def poll(self):
		if self.isComplete() is True:
			return 0
		ret = self.process.poll()
		if ret is not None:
			self.setStatus(ProgressMonitor.StatusCodes.COMPLETED)
		return ret
		
	def communicate(self):
		if self.isComplete() is True:
			return 0
		(stdoutdata, stderrdata) = self.process.communicate()
		ret = self.process.returncode
		self.setStatus(ProgressMonitor.StatusCodes.COMPLETED)
		return ret
#end class ProcessWithProgress			
	
class SmartWorkQueue:
	q = None
	pm = None
	threads = []
	
	def __init__(self, monitor, num_workers):
		self.pm = monitor
		self.q = Queue(0) #no maxsize
		for _ in range(num_workers):
			self.threads.append(SmartWorkQueueWorker(self.q, self.pm))
		
	def addTask(self, cmd, rollback, **kwargs):
		#print "about to add task"
		self.q.put_nowait((cmd, rollback, kwargs))
		#print "Queue size now: "+str(self.size())
	
	def size(self):
		return self.q.qsize()
	
	def empty(self):
		return self.q.empty()
	
	def waitCompletion(self):
		"""Wait for completion of all the tasks in the queue"""
		self.q.join()

	def destruct(self):
		for t in self.threads:
			t.shutdown()
			t.join()
#end class SmartWorkQueue		
		
class SmartWorkQueueWorker(Thread):
	"""Thread executing tasks from a given tasks queue"""
	def __init__(self, tasks, monitor):
		Thread.__init__(self)
		self.shouldExit = Event()
		self.shouldRetry = True
		self.tasks = tasks
		self.monitor = monitor
		self.daemon = True
		self.start()

	def shutdown(self):
		self.shouldExit.set()

	def isShutdown(self):
		return self.shouldExit.isSet()

	def run(self):
		self.monitor = self.monitor.factory()
		while not self.isShutdown():
			try:	
				cmd, rollback, kwargs = self.tasks.get(True, 1)
				p = self.monitor.PopenWithProgress(cmd, rollback, **kwargs)
				p.communicate()
			except Empty:
				continue
			except Exception, e:
				print e
				if self.shouldRetry is True:
					print "re-enqueueing....."
					self.tasks.put_nowait((cmd, rollback, kwargs))
				self.tasks.task_done()
			else:
				self.tasks.task_done()	
#end class SmartWorkQueueWorker(Thread)


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

from tabulate import tabulate
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
		return end - start

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
		return end - start

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

















	
