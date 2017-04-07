import ast
import os
from ThackTech.Pipelines import PipelineModule
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class ReadOutputManifest(PipelineModule):
	
	def __init__(self, **kwargs):
		super(ReadOutputManifest, self).__init__('ReadOutputManifest', 'Read Output Manifest', **kwargs)
	#end __init__()
	
	def run(self, cxt):
		cxt.log.write("-> Reading cxt.sample output manifest....\n")
		cxt.log.flush()
		count = 0
		man_files = {}
		with open(os.path.join(cxt.sample.dest, cxt.sample.name+'_output_manifest.tsv'), 'r') as output_manifest:
			for line in output_manifest:
				line = line.strip()
				if line == '':
					continue
				if line.startswith('id'):
					continue #header line
				parts = line.split('\t')
				if len(parts) > 0:
					c = FileContext(parts[2], int(parts[3]), parts[4], parts[5])
					f = FileInfo(parts[6], c)
					if len(parts) > 7 and parts[7].strip() != '':
						tuples = [item.split("=") for item in parts[7].split(";")]
						for t in tuples:
							f.attributes[t[0]] = ast.literal_eval(t[1])
					if parts[1] == 'None':
						#this is a primary file
						man_files[parts[0]] = f
					else:
						#this is a companion file
						man_files[parts[1]].companions.append(f)

					count += 1
		for f in man_files.values():
			cxt.sample.add_file(f)
		cxt.log.write("-> Read %d items from output manifest....\n" % (count,))
		cxt.log.flush()
	#end run()
#end class ReadOutputManifest