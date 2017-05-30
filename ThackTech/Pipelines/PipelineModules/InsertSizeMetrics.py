import os
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule
import subprocess
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class InsertSizeMetrics(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='InsertSizeMetrics', short_description='Calculating Insert Size Metrics')
		super_args.update(**kwargs)
		super(InsertSizeMetrics, self).__init__(**super_args)
		
		self._name_resolver('bam')
	#end __init__()
	
	def tool_versions(self):
		return {
			"picard-tools:CollectInsertSizeMetrics": self._call_output(["picard-tools", "CollectInsertSizeMetrics", "--version"], stderr=subprocess.STDOUT)
		}
	
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
			insmet_outbase = os.path.join(ism_dir, bam.basename_with_ext('insertsize'))
			args = [
				'picard-tools',
				'CollectInsertSizeMetrics',
				'INPUT={}'.format(bam.fullpath),
				'OUTPUT={}'.format(insmet_outbase+'.txt'),
				'HISTOGRAM_FILE={}'.format(insmet_outbase+'.hist.pdf')
			]
			self._run_subprocess(args, stderr=subprocess.STDOUT, stdout=cxt.log)
			cxt.log.write('\t-> Completed insert size analysis...\n')
			cxt.log.flush()
		
		out_files = []
		out_files.append(FileInfo(insmet_outbase+'.hist.pdf', FileContext.from_module_context(cxt, 'histogram_figure')))
		out_files.append(FileInfo(insmet_outbase+'.txt', FileContext.from_module_context(cxt, 'histogram_data')))
		return out_files
	#end run()
#end InsertSizeMetrics