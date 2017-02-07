import os
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class InsertSizeMetrics(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'InsertSizeMetrics', 'Calculating Insert Size Metrics')
	#end __init__()

	def supported_types(self):
		return ['bam']
	#end supported_types()
	
	def run(self, sample, logfile):
		if not sample.get_attribute('PE'):
			logfile.write('\t-> Skipping insert size analysis because data is not paired-end...\n')
		else:
			logfile.write('\t-> Running insert size analysis...\n')
			logfile.flush()
			ism_dir = os.path.join(sample.dest, 'ism')
			Common.ensure_dir(ism_dir)
			insmet_outbase = os.path.join(ism_dir, sample.name+'.insertsize')
			self._run_subprocess('picard-tools CollectInsertSizeMetrics INPUT=%s OUTPUT=%s HISTOGRAM_FILE=%s' % (bam, insmet_outbase+'.txt', insmet_outbase+'.hist.pdf'), shell=True)
			logfile.write('\t-> Completed insert size analysis...\n')
			logfile.flush()
		
		return {
			'histogram_figure': insmet_outbase+'.hist.pdf',
			'histogram_data':	insmet_outbase+'.txt'
		}
	#end run()
#end InsertSizeMetrics