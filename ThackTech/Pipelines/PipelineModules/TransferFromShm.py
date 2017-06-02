import os
import shutil
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule


class TransferFromShm(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='fromshm', short_description='Transfer From RAM disk')
		super_args.update(**kwargs)
		super(TransferFromShm, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		pass
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		pass
	#end __declare_resolvers()
	
	def run(self, cxt):
		cxt.log.write('-> Moving output to final destination...\n')
		cxt.log.flush()
		
		#first remove the input files, so we can do an easy move....
		for fileinfo in cxt.sample.files:
			if 'origional_location' in fileinfo.attributes:
				os.remove(fileinfo.fullpath)
				fileinfo._set_path(fileinfo.attributes['origional_location'])
				
			for companion in fileinfo.companions:
				if 'origional_location' in companion.attributes:
					os.remove(companion.fullpath)
					companion._set_path(companion.attributes['origional_location'])
		
		
		#first remove the input files, so we can do an easy move....
# 		for key, value in cxt.sample.get_file_group('source').iteritems():
# 			if value.endswith('.bam'):
# 				#dont forget about the index file!
# 				if os.path.isfile(value+'.bai'):
# 					os.remove(value+'.bai') 
# 				elif os.path.isfile(os.path.splitext(value)[0]+'.bai'):
# 					os.remove(os.path.splitext(value)[0]+'.bai')
# 			os.remove(value)
			

		true_dest = cxt.sample.get_attribute('origional_dest')
		filetools.ensure_dir(true_dest)
		procs = []
		

		self._run_subprocess('cp -pr -t '+true_dest+' '+os.path.join(cxt.sample.dest, '*'), shell=True, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		for fileinfo in cxt.sample.files:
			if fileinfo.dirname == 
		
		
# 		for fileinfo in cxt.sample.files:
# 			if groupkey == 'source':
# 				orig_sources = cxt.sample.get_attribute('origional_sources')
# 				for label in orig_sources:
# 					cxt.sample.add_file('source', label, orig_sources[label])
# 			else:
# 				for labelkey in cxt.sample.files[groupkey].keys():
# 					cxt.sample.add_file(groupkey, labelkey, os.path.join(true_dest, os.path.basename(cxt.sample.files[groupkey][labelkey])))
# 		
		shutil.rmtree(cxt.sample.dest) #cleanup shm
		cxt.sample.dest = true_dest
		cxt.sample.remove_attribute('origional_dest')
	#end run()
#end class TransferToShm