import os
import shutil
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule


class TransferFromShm(PipelineModule):
	
	def __init__(self):
		super(TransferFromShm, self).__init__('fromshm', 'Transfer From RAM disk')
	#end __init__()
	
	def run(self, cxt):
		cxt.log.write('-> Moving output to final destination...\n')
		cxt.log.flush()
		#first remove the input files, so we can do an easy move....
		for key, value in cxt.sample.get_file_group('source').iteritems():
			if value.endswith('.bam'):
				#dont forget about the index file!
				if os.path.isfile(value+'.bai'):
					os.remove(value+'.bai') 
				elif os.path.isfile(os.path.splitext(value)[0]+'.bai'):
					os.remove(os.path.splitext(value)[0]+'.bai')
			os.remove(value)
			

		true_dest = cxt.sample.get_attribute('origional_dest')
		filetools.ensure_dir(true_dest)
		self._run_subprocess('cp -pr -t '+true_dest+' '+os.path.join(cxt.sample.dest, '*'), shell=True, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		for groupkey in cxt.sample.files.keys():
			if groupkey == 'source':
				orig_sources = cxt.sample.get_attribute('origional_sources')
				for label in orig_sources:
					cxt.sample.add_file('source', label, orig_sources[label])
			else:
				for labelkey in cxt.sample.files[groupkey].keys():
					cxt.sample.add_file(groupkey, labelkey, os.path.join(true_dest, os.path.basename(cxt.sample.files[groupkey][labelkey])))
		
		shutil.rmtree(cxt.sample.dest) #cleanup shm
		cxt.sample.dest = true_dest
	#end run()
#end class TransferToShm