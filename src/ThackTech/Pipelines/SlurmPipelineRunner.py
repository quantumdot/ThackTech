import os
import sys
import time
import traceback
#from multiprocessing import Pool, Manager 
#from multiprocessing import Manager
import dill	#use dill for pickling, actually supports serializing useful things! (i.e. lambdas, objects)
import multiprocess as mp	#use this ls alternative multiprocessing from pathos, used in combination with dill
import subprocess
from ThackTech.Pipelines import PipelineRunner, GLOBAL_MANAGER
from ThackTech.Processes import MultiStatusProgressItem, MultiStatusProgressBar

try:
	cpu_count = mp.cpu_count()
except:
	cpu_count = 1
if cpu_count is None:
	cpu_count = 1


class SlurmPipelineRunner(PipelineRunner):
	def __init__(self, pipeline, partition="main", nodes=1, threads_per_node=cpu_count, time_limit="1:00:00"):
		PipelineRunner.__init__(self, pipeline)
		self.partition = partition
		self.nodes = nodes
		self.threads_per_node = threads_per_node
		self.time_limit = time_limit

	def run(self, samples):
		curr_time = int(time.time())
		
		with open(os.path.abspath("%s_%d.log" % (self.pipeline.safe_name(), curr_time)), 'w', 0) as logout:
			sample_count = len(samples)
			if sample_count < 1:
				logout.write("No samples to run!")
				return #No samples to run!

			pipeline_pickles = os.path.abspath("pipeline_%s_%d.dill" % (self.pipeline.safe_name(), curr_time))
			with open(pipeline_pickles, 'wb') as f:
				dill.dump(self.pipeline, f)
			
			sample_pickles = []
			status_pickles = []
			self.tasks_statuses = GLOBAL_MANAGER.dict()
			for i in range(len(samples)):
				sample_pickles.append(os.path.abspath("pipeline_%s_%d_s%d.dill" % (self.pipeline.safe_name(), curr_time, i)))
				status_pickles.append(os.path.abspath("pipeline_%s_%d_s%d_status.dill" % (self.pipeline.safe_name(), curr_time, i)))
				with open(sample_pickles[i], 'wb') as sf:
					dill.dump(samples[i], sf)
				with open(status_pickles[i], 'wb') as ssf:
					self.tasks_statuses[samples[i].name] = MultiStatusProgressItem(samples[i].name, 'Queued...', order=i)
					dill.dump(self.tasks_statuses[samples[i].name], ssf)
			try:
				srun_base = [
					'srun',
					'--partition', self.partition,
					'-n', str(self.nodes),
					'--cpus-per-task', str(self.threads_per_node),
					'--time', self.time_limit,
					#'--export', 'ALL',
					'python', os.path.join(os.path.dirname(os.path.abspath(__file__)), "PipelineEntry.py"),
					pipeline_pickles
				]
				for i in range(len(samples)):
					subprocess.Popen(srun_base + [sample_pickles[i], status_pickles[i]], stderr=subprocess.STDOUT, stdout=logout)
				
				progress = MultiStatusProgressBar(sample_count, "Total Progress", barlength=50, handle=sys.stderr).start()
				progress.update(0, None, self.tasks_statuses)
				
				results = []
				while len(results) < sample_count:
					time.sleep(0.5)
					for i in range(len(samples)):
						try:
							with open(status_pickles[i], 'rb') as f:
								self.tasks_statuses[samples[i].name] = dill.load(f)
							if self.tasks_statuses[samples[i].name].is_finished() and i not in results:
								logout.write("task %d is complete!\n" % (i,))
								logout.flush()
								results.append(i)
						except:# EOFError:
							# It is very likely that we get some file IO race conditions here....
							# we dont really care if we get an error too much, because we can probably grab the status
							# on the next pass around. Therefore, just ignore the error.
							#traceback.print_exc(file=logout)
							pass
					progress.update(len(results), '(%d/%d complete)' % (len(results), sample_count), self.tasks_statuses)
				
				progress.update(sample_count, 'Done!', self.tasks_statuses)
				progress.finish()
			except:
				traceback.print_exc(file=logout)
			finally:
				for p in [pipeline_pickles]+sample_pickles+status_pickles:
					if os.path.exists(p):
						subprocess.call(['rm', '-rf', p], stderr=subprocess.STDOUT, stdout=logout)
#end class SlurmPipelineRunner