import os
import sys
from Queue import Queue, Empty
from ThackTech.Processes import SmartWorkQueueWorker


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