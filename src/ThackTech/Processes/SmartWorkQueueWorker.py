import os
import sys
from Queue import Empty
from threading import Thread, Event



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