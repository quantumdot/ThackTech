import sys
from subprocess import PIPE
from ThackTech.Pipelines import PipelineModule


class DumpEnv(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='DumpEnv', short_description='Dump environment information')
		super_args.update(**kwargs)
		super(DumpEnv, self).__init__(**super_args)
	#end __init__()
	
	def __declare_parameters(self):
		pass
		#self.add_parameter(ModuleParameter('sleep_time', int, 10, desc="Amount of time, in seconds, to sleep after saying hello."))
	#end __declare_parameters()
	
	def __declare_resolvers(self):
		pass
	#end __declare_resolvers()
	
	def load_modules(self, log):
		# """ Loads any system modules required for this module to function. """
		# log.write("loading samtools")
		# o, e = self._run_subprocess("module load samtools/0.1.19", shell=True, stdout=PIPE, stderr=PIPE)
		# log.write(str(o)+"\n"+str(e))
		# log.flush()
		
		try:
			#from ThackTech.Pipelines import lmodHelper
			#log.write(lmodHelper.module("load", "samtools/0.1.19"))
			
			
# 			from env_modules_python import module
# 			log.write(module("load", "samtools/0.1.19"))
# 			log.flush()

			from ThackTech.Pipelines import LmodHelper2
			log.write(LmodHelper2.module("load", "samtools/0.1.19"))
			log.flush()
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
		
		cxt.log.write("Getting python version information...\n")
		cxt.log.write(sys.version)
		
		self.load_modules(cxt.log)
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