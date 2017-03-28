import os
from ThackTech.Pipelines import PipelineModule


class GeneratePlainBed(PipelineModule):
	def __init__(self, **kwargs):
		super(GeneratePlainBed, self).__init__('PlainBed', 'Generate bed from peaks', **kwargs)
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, cxt):
		try_paths = [
			os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.narrowPeak'),
			os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.broadPeak'),
			os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.gappedPeak')
		]
		
		for path in try_paths:
			if os.path.isfile(path):
				with open(path, 'r') as encodepeaks:
					with open(os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.bed'), 'w') as outfile:
						self._run_subprocess(['cut', '-f', '1-6'], stdin=encodepeaks, stderr=cxt.log, stdout=outfile, cwd=cxt.sample.dest)
				return {
					'simple_bed': os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.bed')
				}
		cxt.log.write('no conversion necessary...\n')
		cxt.log.flush()
	#end run()
#end class GeneratePlainBed