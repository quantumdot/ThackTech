from ThackTech.Pipelines.PipelineRunner import PipelineRunner
from ThackTech.Pipelines.PipelineRunner import _execute_pipeline_on_sample
from ThackTech.Processes import MultiStatusProgressItem

class SerialPipelineRunner(PipelineRunner):
	def __init__(self, pipeline):
		PipelineRunner.__init__(self, pipeline)

	def run(self, samples):
		for sample in samples:
			_execute_pipeline_on_sample(self.pipeline, sample, {sample.name: MultiStatusProgressItem(sample.name, 'Queued...')})
#end class SerialPipelineRunner