import filecmp
import os
import shlex
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class TransferToShm(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'toshm', 'Transfer To RAM disk')
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, sample, logfile):
		logfile.write('-> Copying source files to RAM disk....\n')
		paths = sample.get_file_group('source')
		sample.add_attribute('origional_sources', paths)
		sample.add_attribute('origional_dest', sample.dest)
		shmdest = os.path.join(self.get_parameter_value('shm_path'), sample.name)
		filetools.ensure_dir(shmdest)
		
		logfile.flush()
		procs = []
		newpaths = {}
		for key in paths:
			fullpath = paths[key]
			filename = os.path.basename(fullpath)
			dir = os.path.dirname(fullpath)
			basename, ext = os.path.splitext(filename)
			logfile.write('-> preparing to move %s\n' % (fullpath,))
			
			dest = os.path.join(shmdest, filename)
			newpaths[key] = dest
			#check if the file exists, and if so do the files look alike
			#if the file exists and looks the same, we can skip copying it over (but we still return the shm path) 
			if os.path.isfile(dest) and filecmp.cmp(fullpath, dest):
				logfile.write('    -> It appears that the file already exists on RAMFS. Skipping copying of this file\n')
				logfile.flush()
			else:
				if ext == '.bam':
					#dont forget the index!
					if os.path.isfile(fullpath+'.bai'):
						procs.append(subprocess.Popen(['cp', '-p', '-t', shmdest, fullpath+'.bai'], stderr=subprocess.STDOUT, stdout=logfile))
					elif os.path.isfile(os.path.join(dir, basename+'.bai')):
						procs.append(subprocess.Popen(['cp', '-p', '-t', shmdest, os.path.join(dir, basename+'.bai')], stderr=subprocess.STDOUT, stdout=logfile))
					
				procs.append(subprocess.Popen(['cp', '-p', '-t', shmdest, fullpath], stderr=subprocess.STDOUT, stdout=logfile))
			
			sample.add_file('source', key, dest)
		if len(procs) > 0:
			for p in procs:
				p.communicate()
		
		sample.dest = shmdest
	#end run()
#end class TransferToShm