import os
import sys
import sqlite3
from ThackTech.Processes import ProgressMonitor, ProgressMonitorAdapter, ProgressStatusCodes




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
			return ProgressStatusCodes.UNKNOWN
		return r[0]
		
	def setStatus(self, task, status):
		c = self.store.cursor()
		c.execute('INSERT OR REPLACE INTO '+self.table_name+' ('+self.task_col+', '+self.stat_col+') VALUES (?, ?)', (task, status))
		self.store.commit()
#end class SQLiteProgressMonitorAdapter