import os
import sys
from ThackTech.Processes import ProgressStatusCodes




class ProcessTicket:
	
	def __init__(self, signature, status=ProgressStatusCodes.UNKNOWN):
		self.times = {}
		self.signature = signature
		self.status = status
#end class ProcessTicket