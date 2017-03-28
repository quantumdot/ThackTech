import time
import dill	#use dill for pickling, actually supports serializing useful things! (i.e. lambdas, objects)
import multiprocess as mp	#use this ls alternative multiprocessing from pathos, used in combination with dill
from ThackTech.Pipelines import PipelineRunner, GLOBAL_MANAGER, CPU_COUNT
from ThackTech.Pipelines.PipelineRunner import _execute_pipeline_on_sample
from ThackTech.Processes import MultiStatusProgressItem, MultiStatusProgressBar


class ParallelPipelineRunner(PipelineRunner):
	"""Run a pipeline against samples on a single node in parallel
	
	This PipelineRunner uses the multiprocess worker pool to parallelize running
	of pipelines on multiple samples.
	"""
	def __init__(self, pipeline, threads=CPU_COUNT):
		"""Initialize this PipelineRunner
		
		Parameters:
			pipeline:	(AnalysisPipeline)	Pipeline to run
			threads:	(int)	Number of workers to spawn
		"""
		PipelineRunner.__init__(self, pipeline)
		self.threads = threads
	#end __init__()

	def run(self, samples):
		"""Run this pipeline runner on the supplied samples
		
		Parameters:
			samples: (iterable of PipelineSample) samples to run against
		"""
		sample_count = len(samples)
		if sample_count < 1:
			return #No samples to run!
		
		self.tasks_statuses = GLOBAL_MANAGER.dict()
		for sname in [s.name for s in samples]:
			self.tasks_statuses[sname] = MultiStatusProgressItem(sname, 'Queued...')
		pool = mp.Pool(processes=self.threads)
		
		progress = MultiStatusProgressBar(sample_count, "Total Progress", barlength=50).start()
		progress.update(0, None, self.tasks_statuses)
		
		results = []
		r = [pool.apply_async(_execute_pipeline_on_sample, (self.pipeline, sample, self.tasks_statuses), callback=results.append) for sample in samples]
		
		while len(results) < sample_count:
			time.sleep(0.5)
			progress.update(len(results), '(%d/%d complete)' % (len(results), sample_count), self.tasks_statuses)
		
		pool.close()
		pool.join()
		progress.update(sample_count, 'Done!', self.tasks_statuses)
		progress.finish()
	#end run()
#end class ParallelPipelineRunner