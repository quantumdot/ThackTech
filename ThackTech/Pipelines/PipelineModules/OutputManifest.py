import os
from ThackTech.Pipelines import PipelineModule


class OutputManifest(PipelineModule):
	
	def __init__(self, **kwargs):
		super(OutputManifest, self).__init__('OutputManifest', 'Generate Output Manifest', **kwargs)
	#end __init__()

	
	def run(self, cxt):
		cxt.log.write("-> Writing cxt.sample output manifest....\n")
		cxt.log.flush()
		dest = os.path.join(cxt.sample.dest, cxt.sample.name+'_output_manifest.tsv')
		if os.path.exists(dest):
			write_headers = False
		else:
			write_headers = True
		
		with open(dest, 'a') as output_manifest:
			if write_headers:
				output_manifest.write('module\tlabel\tpath\n')
			for groupkey in cxt.sample.files.keys():
				for labelkey in cxt.sample.files[groupkey].keys():
					val = cxt.sample.files[groupkey][labelkey]
					if isinstance(val, str):
						out_val = "'%s'" % (val,)
					else:
						out_val = str(val)
					output_manifest.write('%s\t%s\t%s\n' % (groupkey, labelkey, out_val))
	#end run()
#end class OutputManifest

