import os
import platform
from subprocess import PIPE
import sys
import time

import ThackTech.Common as Common
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class DumpEnv(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'DumpEnv', 'Dump environment information')
		
		#self.add_parameter(ModuleParameter('sleep_time', int, 10, desc="Amount of time, in seconds, to sleep after saying hello."))
	#end __init__()
	
	def load_modules(self, logfile):
		# """ Loads any system modules required for this module to function. """
		# logfile.write("loading samtools")
		# o, e = self._run_subprocess("module load samtools/0.1.19", shell=True, stdout=PIPE, stderr=PIPE)
		# logfile.write(str(o)+"\n"+str(e))
		# logfile.flush()
		
		try:
			#from ThackTech.Pipelines import lmodHelper
			#logfile.write(lmodHelper.module("load", "samtools/0.1.19"))
			from env_modules_python import module
			logfile.write(module("load", "samtools/0.1.19"))
			logfile.flush()
		except:
			pass
	#end load_modules()
	
	def run(self, sample, logfile):
		# logfile.write("Dumping Python Environment Information (os.environ):\n")
		# logfile.write("%s\n\n\n" % (os.environ,))
		# logfile.flush()
		# logfile.write("Dumping bash Environment Information (bash env):\n")
		# o, e = self._run_subprocess(['env'], stdout=PIPE, stderr=PIPE)
		# logfile.write(str(o)+"\n"+str(e))
		# logfile.flush()
		
		# try:
			# #from ThackTech.Pipelines import lmodHelper
			# #logfile.write(lmodHelper.module("load", "samtools/0.1.19"))
			# from env_modules_python import module
			# logfile.write(module("load", "samtools/0.1.19"))
			# logfile.flush()
		# except:
			# pass
			
		# logfile.write("Dumping Python Environment Information (os.environ):\n")
		# logfile.write("%s\n\n\n" % (os.environ,))
		# logfile.flush()
		# logfile.write("Dumping bash Environment Information (bash env):\n")
		# o, e = self._run_subprocess(['env'], stdout=PIPE, stderr=PIPE)
		# logfile.write(str(o)+"\n"+str(e))
		# logfile.flush()
		
		#logfile.write("Dumping inline Environment Information \n")
		#logfile.flush()
		logfile.write("running samtools \n")
		logfile.flush()
		o, e = self._run_subprocess('module list', shell=True, stdout=PIPE, stderr=PIPE)
		logfile.write(str(o)+"\n"+str(e))
		logfile.flush()
		#self._run_subprocess(['samtools'])
		logfile.flush()
		#self._run_subprocess('echo "current environment\n"; env; echo "module loading\n"; module load samtools/0.1.19; echo "new environment\n"; env;', shell=True, stdout=PIPE, stderr=PIPE, executable='/bin/bash')
		#logfile.write(str(o)+"\n"+str(e))
		#logfile.flush()
	#end run()
#end class DumpEnv