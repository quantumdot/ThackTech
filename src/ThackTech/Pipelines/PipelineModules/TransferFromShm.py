import os
import shutil
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class TransferFromShm(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'fromshm', 'Transfer From RAM disk')
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, sample, logfile):
		logfile.write('-> Moving output to final destination...\n')
		logfile.flush()
		#first remove the input files, so we can do an easy move....
		for key, value in sample.get_file_group('source').iteritems():
			if value.endswith('.bam'):
				#dont forget about the index file!
				if os.path.isfile(value+'.bai'):
					os.remove(value+'.bai') 
				elif os.path.isfile(os.path.splitext(value)[0]+'.bai'):
					os.remove(os.path.splitext(value)[0]+'.bai')
			os.remove(value)
			

		true_dest = sample.get_attribute('origional_dest')
		Common.ensure_dir(true_dest)
		self._run_subprocess('cp -pr -t '+true_dest+' '+os.path.join(sample.dest, '*'), shell=True, stderr=subprocess.STDOUT, stdout=logfile)
		
		for groupkey in sample.files.keys():
			if groupkey == 'source':
				orig_sources = sample.get_attribute('origional_sources')
				for label in orig_sources:
					sample.add_file('source', label, orig_sources[label])
			else:
				for labelkey in sample.files[groupkey].keys():
					sample.add_file(groupkey, labelkey, os.path.join(true_dest, os.path.basename(sample.files[groupkey][labelkey])))
		
		shutil.rmtree(sample.dest) #cleanup shm
		sample.dest = true_dest
	#end run()
#end class TransferToShm