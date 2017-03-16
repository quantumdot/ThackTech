import ast
import filecmp
import os
import shlex
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class ReadOutputManifest(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'ReadOutputManifest', 'Read Output Manifest')
	#end __init__()
	
	def run(self, sample, logfile):
		logfile.write("-> Reading sample output manifest....\n")
		logfile.flush()
		count = 0
		with open(os.path.join(sample.dest, sample.name+'_output_manifest.tsv'), 'r') as output_manifest:
			for line in output_manifest:
				line = line.strip()
				if line == '':
					continue
				if line == 'module\tlabel\tpath':
					continue
				parts = line.split('\t')
				if len(parts) > 0:
					sample.add_file(parts[0], parts[1], ast.literal_eval(parts[2]))
					count += 1
		logfile.write("-> Read %d items from output manifest....\n" % (count,))
		logfile.flush()
	#end run()
#end class ReadOutputManifest