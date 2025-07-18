import os
import sys
import time
import traceback
#from multiprocessing import Pool, Manager 
#from multiprocessing import Manager
import dill	#use dill for pickling, actually supports serializing useful things! (i.e. lambdas, objects)
import multiprocess as mp	#use this alternative multiprocessing from pathos, used in combination with dill
import subprocess
import uuid
from ThackTech import filetools
from ThackTech.Pipelines import PipelineRunner, GLOBAL_MANAGER, CPU_COUNT
from ThackTech.Processes import MultiStatusProgressItem, MultiStatusProgressBar


class SlurmPipelineRunner(PipelineRunner):
	"""Run a pipeline against samples on a SLURM managed cluster/system
	
	This PipelineRunner uses the multiprocess worker pool to parallelize running
	of pipelines on multiple samples.
	"""
	def __init__(self, pipeline, partition="main", nodes=1, threads=CPU_COUNT, time_limit="1:00:00", mem='4G', slurm_cmd='srun'):
		"""
		
		Parameters:
			pipeline:	(AnalysisPipeline)	Pipeline to run
			partition:	(string)			System partition to request for execution
			nodes:		(int) 				Number of nodes to request per sample
			threads:	(int)				Number of CPUs to request per sample
			time_limit:	(time-string)		Approximate time limit for this SLURM job to run.
			slurm_cmd:	(string)			Slurm command to use when spawning jobs, one of [srun, sbatch]
		"""
		PipelineRunner.__init__(self, pipeline)
		self.slurm_cmd = slurm_cmd
		self.partition = partition
		self.nodes = nodes
		self.threads_per_node = threads
		self.time_limit = time_limit
		self.mem = mem

	def run(self, samples):
		curr_time = int(time.time())
		uid = uuid.uuid4()
		
		run_dir = os.path.abspath(".pj/")
		filetools.ensure_dir(run_dir)
		prefix = "{pipename}_{uid}".format(pipename=self.pipeline.safe_name, uid=curr_time)

		with open(os.path.join(run_dir, "{prefix}.log".format(prefix=prefix)), 'w', 0) as logout:
			sample_count = len(samples)
			if sample_count < 1:
				logout.write("No samples to run!")
				logout.flush()
				return #No samples to run!

			pipeline_pickles = os.path.join(run_dir, "{prefix}.dill".format(prefix=prefix))
			with open(pipeline_pickles, 'wb') as f:
				dill.dump(self.pipeline, f)
			
			sample_pickles = []
			status_pickles = []
			self.tasks_statuses = GLOBAL_MANAGER.dict()
			for i in range(len(samples)):
				sample_prefix = "{prefix}_s{index}_{sname}".format(prefix=prefix, index=i, sname=samples[i].name)
				
				sample_pickles.append(os.path.join(run_dir, sample_prefix+".dill"))
				status_pickles.append(os.path.join(run_dir, sample_prefix+"_status.dill"))
				with open(sample_pickles[i], 'wb') as sf:
					dill.dump(samples[i], sf)
				with open(status_pickles[i], 'wb') as ssf:
					self.tasks_statuses[samples[i].name] = MultiStatusProgressItem(samples[i].name, 'Queued...', order=i)
					dill.dump(self.tasks_statuses[samples[i].name], ssf)
			try:
				for i in range(len(samples)):
					if self.slurm_cmd == 'srun':
						srun_cmd = [
							'srun',
							'--job-name', sample_prefix,
							'--partition', str(self.partition),
							'-n', str(self.nodes),
							'--cpus-per-task', str(self.threads_per_node),
							'--time', str(self.time_limit),
							'--mem', self.mem,
							#'--export', 'ALL',
							'python', os.path.join(os.path.dirname(os.path.abspath(__file__)), "PipelineEntry.py"),
							pipeline_pickles,
							sample_pickles[i],
							status_pickles[i]
						]
						logout.write("Running srun command:\n")
						logout.write(" ".join(srun_cmd))
						logout.write("\n\n")
						logout.flush()
					else:
						srun_cmd = [
							'sbatch',
							'--job-name', sample_prefix,
							'--partition', str(self.partition),
							'-n', str(self.nodes),
							'--cpus-per-task', str(self.threads_per_node),
							'--time', str(self.time_limit),
							'--mem', self.mem,
							#'--export', 'ALL',
							'--output', os.path.join(run_dir, sample_prefix+".out"),
							'--error', os.path.join(run_dir, sample_prefix+".err"),
							'--wrap', ' '.join([
								'python', os.path.join(os.path.dirname(os.path.abspath(__file__)), "PipelineEntry.py"),
								pipeline_pickles,
								sample_pickles[i],
								status_pickles[i]
							])
						]
						logout.write("Running sbatch command:\n")
						logout.write(" ".join(srun_cmd))
						logout.write("\n\n")
						logout.flush()
					subprocess.Popen(srun_cmd, stderr=subprocess.STDOUT, stdout=logout)
				
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
				logout.flush()
			finally:
				for p in [pipeline_pickles]+sample_pickles+status_pickles:
					if os.path.exists(p):
						subprocess.call(['rm', '-rf', p], stderr=subprocess.STDOUT, stdout=logout)
						logout.flush()
	#end run()
	
	def monitor(self):
		PipelineRunner.monitor(self)					
	#end monitor()
#end class SlurmPipelineRunner