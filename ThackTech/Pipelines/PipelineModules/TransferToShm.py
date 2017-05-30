import filecmp
import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule


class TransferToShm(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='toshm', short_description='Transfer To RAM disk')
		super_args.update(**kwargs)
		super(TransferToShm, self).__init__(**super_args)
		
		self._name_resolver('files')
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, cxt):
		cxt.log.write('-> Copying source files to RAM disk....\n')
		paths = self.resolve_input('files', cxt)
		
		shmdest = os.path.join(self.get_parameter_value('shm_path'), cxt.sample.name)
		filetools.ensure_dir(shmdest)
		cxt.sample.set_attribute('origional_dest', cxt.sample.dest)
		cxt.sample.dest = shmdest
		
		cxt.log.flush()
		procs = []
		for fileinfo in paths:			
			cxt.log.write('-> preparing to move %s\n' % (fileinfo.fullpath,))

			#check if the file exists, and if so do the files look alike
			#if the file exists and looks the same, we can skip copying it over (but we still return the shm path) 
			if os.path.isfile(os.path.join(shmdest, fileinfo.basename)) and filecmp.cmp(fileinfo.fullpath, os.path.join(shmdest, fileinfo.basename)):
				cxt.log.write('    -> It appears that the file already exists on RAMFS. Skipping copying of this file\n')
				cxt.log.flush()
			else:
				for companion in fileinfo.companions:
					procs.append(subprocess.Popen(['cp', '-p', '-t', shmdest, companion.fullpath], stderr=subprocess.STDOUT, stdout=cxt.log))
					companion.attributes['origional_location'] = companion.fullpath
					companion._set_path(os.path.join(shmdest, companion.basename))
					
				procs.append(subprocess.Popen(['cp', '-p', '-t', shmdest, fileinfo.fullpath], stderr=subprocess.STDOUT, stdout=cxt.log))
				fileinfo.attributes['origional_location'] = fileinfo.fullpath
				fileinfo._set_path(os.path.join(shmdest, fileinfo.basename))
			
		if len(procs) > 0:
			for p in procs:
				p.communicate()
		
	#end run()
#end class TransferToShm