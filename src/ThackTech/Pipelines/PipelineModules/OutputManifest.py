import filecmp
import os
import shlex
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class OutputManifest(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'OutputManifest', 'Generate Output Manifest')
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, sample, logfile):
		logfile.write("-> Writing sample output manifest....\n")
		logfile.flush()
		dest = os.path.join(sample.dest, sample.name+'_output_manifest.tsv')
		if os.path.exists(dest):
			write_headers = False
		else:
			write_headers = True
		
		with open(dest, 'a') as output_manifest:
			if write_headers:
				output_manifest.write('module\tlabel\tpath\n')
			for groupkey in sample.files.keys():
				for labelkey in sample.files[groupkey].keys():
					val = sample.files[groupkey][labelkey]
					if isinstance(val, str):
						out_val = "'%s'" % (val,)
					else:
						out_val = str(val)
					output_manifest.write('%s\t%s\t%s\n' % (groupkey, labelkey, out_val))
	#end run()
#end class OutputManifest

