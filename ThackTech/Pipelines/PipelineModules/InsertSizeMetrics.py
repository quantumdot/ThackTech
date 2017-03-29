import os
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule


class InsertSizeMetrics(PipelineModule):
	
	def __init__(self, **kwargs):
		super(InsertSizeMetrics, self).__init__('InsertSizeMetrics', 'Calculating Insert Size Metrics', **kwargs)
		
		self._name_resolver('bam')
	#end __init__()
	
	def run(self, cxt):
		if not cxt.sample.get_attribute('PE'):
			cxt.log.write('\t-> Skipping insert size analysis because data is not paired-end...\n')
			return
		else:
			cxt.log.write('\t-> Running insert size analysis...\n')
			cxt.log.flush()
			bam = self.resolve_input('bam', cxt)
			ism_dir = os.path.join(cxt.sample.dest, 'ism')
			filetools.ensure_dir(ism_dir)
			insmet_outbase = os.path.join(ism_dir, cxt.sample.name+'.insertsize')
			args = [
				'picard-tools',
				'CollectInsertSizeMetrics',
				'INPUT={}'.format(bam),
				'OUTPUT={}'.format(insmet_outbase+'.txt'),
				'HISTOGRAM_FILE={}'.format(insmet_outbase+'.hist.pdf')
			]
			self._run_subprocess(args)
			cxt.log.write('\t-> Completed insert size analysis...\n')
			cxt.log.flush()
		
		return {
			'histogram_figure': insmet_outbase+'.hist.pdf',
			'histogram_data':	insmet_outbase+'.txt'
		}
	#end run()
#end InsertSizeMetrics