import ast
import os
from ThackTech.Pipelines import PipelineModule


class ReadOutputManifest(PipelineModule):
	
	def __init__(self, **kwargs):
		super(ReadOutputManifest, self).__init__('ReadOutputManifest', 'Read Output Manifest', **kwargs)
	#end __init__()
	
	def run(self, cxt):
		cxt.log.write("-> Reading cxt.sample output manifest....\n")
		cxt.log.flush()
		count = 0
		with open(os.path.join(cxt.sample.dest, cxt.sample.name+'_output_manifest.tsv'), 'r') as output_manifest:
			for line in output_manifest:
				line = line.strip()
				if line == '':
					continue
				if line == 'module\tlabel\tpath':
					continue
				parts = line.split('\t')
				if len(parts) > 0:
					cxt.sample.add_file(parts[0], parts[1], ast.literal_eval(parts[2]))
					count += 1
		cxt.log.write("-> Read %d items from output manifest....\n" % (count,))
		cxt.log.flush()
	#end run()
#end class ReadOutputManifest