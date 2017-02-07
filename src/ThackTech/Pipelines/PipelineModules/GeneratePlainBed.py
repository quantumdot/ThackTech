import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class GeneratePlainBed(PipelineModule):
	def __init__(self):
		PipelineModule.__init__(self, 'PlainBed', 'Generate bed from peaks')
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, sample, logfile):
		try_paths = [
			os.path.join(sample.dest, sample.name+'_peaks.narrowPeak'),
			os.path.join(sample.dest, sample.name+'_peaks.broadPeak'),
			os.path.join(sample.dest, sample.name+'_peaks.gappedPeak')
		]
		
		for path in try_paths:
			if os.path.isfile(path):
				with open(path, 'r') as encodepeaks:
					with open(os.path.join(sample.dest, sample.name+'_peaks.bed'), 'w') as outfile:
						self._run_subprocess(['cut', '-f', '1-6'], stdin=encodepeaks, stderr=logfile, stdout=outfile, cwd=sample.dest)
				return {
					'simple_bed': os.path.join(sample.dest, sample.name+'_peaks.bed')
				}
		logfile.write('no conversion necessary...\n')
		logfile.flush()
	#end run()
#end class GeneratePlainBed