from subprocess import PIPE
from ThackTech.Pipelines import PipelineModule


class DumpEnv(PipelineModule):
	
	def __init__(self, **kwargs):
		super(DumpEnv, self).__init__('DumpEnv', 'Dump environment information', **kwargs)
		
		#self.add_parameter(ModuleParameter('sleep_time', int, 10, desc="Amount of time, in seconds, to sleep after saying hello."))
	#end __init__()
	
	def load_modules(self, cxt):
		# """ Loads any system modules required for this module to function. """
		# cxt.log.write("loading samtools")
		# o, e = self._run_subprocess("module load samtools/0.1.19", shell=True, stdout=PIPE, stderr=PIPE)
		# cxt.log.write(str(o)+"\n"+str(e))
		# cxt.log.flush()
		
		try:
			#from ThackTech.Pipelines import lmodHelper
			#cxt.log.write(lmodHelper.module("load", "samtools/0.1.19"))
			
			
# 			from env_modules_python import module
# 			cxt.log.write(module("load", "samtools/0.1.19"))
# 			cxt.log.flush()

			from ThackTech.Pipelines import LmodHelper2
			cxt.log.write(LmodHelper2.module("load", "samtools/0.1.19"))
			cxt.log.flush()
		except:
			pass
	#end load_modules()
	
	def run(self, cxt):
		# cxt.log.write("Dumping Python Environment Information (os.environ):\n")
		# cxt.log.write("%s\n\n\n" % (os.environ,))
		# cxt.log.flush()
		# cxt.log.write("Dumping bash Environment Information (bash env):\n")
		# o, e = self._run_subprocess(['env'], stdout=PIPE, stderr=PIPE)
		# cxt.log.write(str(o)+"\n"+str(e))
		# cxt.log.flush()
		
		# try:
			# #from ThackTech.Pipelines import lmodHelper
			# #cxt.log.write(lmodHelper.module("load", "samtools/0.1.19"))
			# from env_modules_python import module
			# cxt.log.write(module("load", "samtools/0.1.19"))
			# cxt.log.flush()
		# except:
			# pass
			
		# cxt.log.write("Dumping Python Environment Information (os.environ):\n")
		# cxt.log.write("%s\n\n\n" % (os.environ,))
		# cxt.log.flush()
		# cxt.log.write("Dumping bash Environment Information (bash env):\n")
		# o, e = self._run_subprocess(['env'], stdout=PIPE, stderr=PIPE)
		# cxt.log.write(str(o)+"\n"+str(e))
		# cxt.log.flush()
		
		#cxt.log.write("Dumping inline Environment Information \n")
		#cxt.log.flush()
		cxt.log.write("running samtools \n")
		cxt.log.flush()
		o, e = self._run_subprocess('module list', shell=True, stdout=PIPE, stderr=PIPE)
		cxt.log.write(str(o)+"\n"+str(e))
		cxt.log.flush()
		#self._run_subprocess(['samtools'])
		cxt.log.flush()
		#self._run_subprocess('echo "current environment\n"; env; echo "module loading\n"; module load samtools/0.1.19; echo "new environment\n"; env;', shell=True, stdout=PIPE, stderr=PIPE, executable='/bin/bash')
		#cxt.log.write(str(o)+"\n"+str(e))
		#cxt.log.flush()
	#end run()
#end class DumpEnv