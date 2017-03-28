import filecmp
import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule


class TransferToShm(PipelineModule):
	
	def __init__(self):
		super(TransferToShm, self).__init__('toshm', 'Transfer To RAM disk')
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, cxt):
		cxt.log.write('-> Copying source files to RAM disk....\n')
		paths = cxt.sample.get_file_group('source')
		cxt.sample.add_attribute('origional_sources', paths)
		cxt.sample.add_attribute('origional_dest', cxt.sample.dest)
		shmdest = os.path.join(self.get_parameter_value('shm_path'), cxt.sample.name)
		filetools.ensure_dir(shmdest)
		
		cxt.log.flush()
		procs = []
		newpaths = {}
		for key in paths:
			fullpath = paths[key]
			filename = os.path.basename(fullpath)
			dir = os.path.dirname(fullpath)
			basename, ext = os.path.splitext(filename)
			cxt.log.write('-> preparing to move %s\n' % (fullpath,))
			
			dest = os.path.join(shmdest, filename)
			newpaths[key] = dest
			#check if the file exists, and if so do the files look alike
			#if the file exists and looks the same, we can skip copying it over (but we still return the shm path) 
			if os.path.isfile(dest) and filecmp.cmp(fullpath, dest):
				cxt.log.write('    -> It appears that the file already exists on RAMFS. Skipping copying of this file\n')
				cxt.log.flush()
			else:
				if ext == '.bam':
					#dont forget the index!
					if os.path.isfile(fullpath+'.bai'):
						procs.append(subprocess.Popen(['cp', '-p', '-t', shmdest, fullpath+'.bai'], stderr=subprocess.STDOUT, stdout=cxt.log))
					elif os.path.isfile(os.path.join(dir, basename+'.bai')):
						procs.append(subprocess.Popen(['cp', '-p', '-t', shmdest, os.path.join(dir, basename+'.bai')], stderr=subprocess.STDOUT, stdout=cxt.log))
					
				procs.append(subprocess.Popen(['cp', '-p', '-t', shmdest, fullpath], stderr=subprocess.STDOUT, stdout=cxt.log))
			
			cxt.sample.add_file('source', key, dest)
		if len(procs) > 0:
			for p in procs:
				p.communicate()
		
		cxt.sample.dest = shmdest
	#end run()
#end class TransferToShm