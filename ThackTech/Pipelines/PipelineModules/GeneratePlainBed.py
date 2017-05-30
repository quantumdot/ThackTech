import os
from ThackTech.Pipelines import PipelineModule


class GeneratePlainBed(PipelineModule):
	def __init__(self, **kwargs):
		super_args = dict(name='PlainBed', short_description='Generate bed6 from peaks')
		super_args.update(**kwargs)
		super(GeneratePlainBed, self).__init__(**super_args)
		
		self._name_resolver('encode_beds')
	#end __init__()
	
	def run(self, cxt):
		try_paths = self.resolve_input('encode_beds', cxt)
		
		ext_to_ranges = {
			'.narrowPeak': '$1, $2, $3, $4, $7, $6',
			'.broadPeak':  '$1, $2, $3, $4, $7, $6',
			'.gappedPeak': '$1, $2, $3, $4, $13, $6'	
		}
		
		output_files = {}
		for path in try_paths:
			outbed = os.path.join(cxt.sample.dest, path.basename_with_ext('bed'))
			
			with open(path.fullpath, 'r') as encodepeaks, open(outbed, 'w') as outfile:
				cmd = ['awk', r'{OFS="\t"; print '+ext_to_ranges[path.ext]+'}']
				cxt.log.write('Converting {} to BED6 format and writing to {}....'.format(path.fullpath, outbed))
				cxt.log.write("\n..............................................\n")
				cxt.log.write(" ".join(cmd))
				cxt.log.write("\n..............................................\n")
				cxt.log.flush()
				self._run_subprocess(cmd, stdin=encodepeaks, stderr=cxt.log, stdout=outfile, cwd=cxt.sample.dest)
				
			output_files[path.cxt.role] = outbed
	#end run()
#end class GeneratePlainBed