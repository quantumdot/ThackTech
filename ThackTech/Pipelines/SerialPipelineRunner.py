from ThackTech.Pipelines.PipelineRunner import PipelineRunner
from ThackTech.Pipelines.PipelineRunner import _execute_pipeline_on_sample
from ThackTech.Processes import MultiStatusProgressItem, MultiStatusProgressBar
import threading
import time

class SerialPipelineRunner(PipelineRunner):
	"""Run a pipeline against samples on a single node serially
	
	This PipelineRunner executes the pipeline on one sample at a time on a single node
	"""
	def __init__(self, pipeline):
		PipelineRunner.__init__(self, pipeline)

	def run(self, samples):
		
		sample_count = len(samples)
		if sample_count < 1:
			return #No samples to run!
		
		self.tasks_statuses = dict()
		for sname in [s.name for s in samples]:
			self.tasks_statuses[sname] = MultiStatusProgressItem(sname, 'Queued...')

		progress = MultiStatusProgressBar(sample_count, "Total Progress", barlength=50).start()
		progress.update(0, None, self.tasks_statuses)
		
		for i in range(sample_count):
			thread = threading.Thread(target=_execute_pipeline_on_sample, args=(self.pipeline, samples[i], self.tasks_statuses))
			thread.start()
		
			while thread.is_alive():
				time.sleep(0.5)
				progress.update(i, '(%d/%d complete)' % (i, sample_count), self.tasks_statuses)
			thread.join()
	#end run()
		
	def monitor(self):
		PipelineRunner.monitor(self)
	#end monitor()
#end class SerialPipelineRunner