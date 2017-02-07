import os
import sys
import traceback
import platform
import time
from ThackTech import filetools


class PipelineRunner(object):
	"""Base class for pipeline runners
	
	Pipeline runners separate the logic of actually executing the pipeline, which may differ
	depending on the environment, from the pipeline definition.
	
	To implement a pipeline runner, you must implement the run() method.
	"""
	def __init__(self, pipeline):
		self.pipeline = pipeline
		self.tasks_statuses = None

	def preload_modules(self):
		for step in self.pipeline.pipeline:
			step.load_modules(sys.stderr)
	#end preload_modules()
	
	def run(self, samples):
		pass
#end class PipelineRunner

def _execute_pipeline_on_sample(pipeline, sample, tasks_statuses):
	tasks_statuses[sample.name] = tasks_statuses[sample.name].start()
	#ensure the output directory exists
	filetools.ensure_dir(sample.dest)

	with open(os.path.join(sample.dest, sample.name+'.log'), 'a', buffering=0) as logfile:
		sys.stdout = sys.stderr = logfile
		logfile.write(pipeline.documentation())
		status_counts = {}
		try:
			logfile.write('Processing sample "%s" [%s]....\n' % (sample.name, sample.format))
			logfile.write('-> Pipeline: %s\n' % (pipeline.name,))
			logfile.write('-> Running on: %s\n' % (' '.join(platform.uname()),))
			logfile.write("-> Wall clock: %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"),))
			logfile.write("--------------------------------------------\n\n")
			tasks_statuses[sample.name] = tasks_statuses[sample.name].update(0, 'Preparing...')
			
			pipeline_size = len(pipeline.pipeline)
			status_counts['total'] = pipeline_size
			status_counts['attempted'] = 0
			status_counts['warn'] = 0
			status_counts['critical'] = 0
			for i in range(pipeline_size):
				step = pipeline.pipeline[i]
				status_counts['attempted'] += 1
				logfile.write('Running pipeline step: %s\n' % (step.name,))
				logfile.write("-> Wall clock: %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"),))
				logfile.flush()
				try:
					if not step.is_compatible_with_sample(sample):
						raise ValueError('%s does not support the "%s" file format!\nThe following types ARE supported: %s\nSkipping processing......\n' % (step.name, sample.format, ', '.join(step.supported_types())))
				
					tasks_statuses[sample.name] = tasks_statuses[sample.name].update(float(i)/float(pipeline_size), step.description)
					#step.load_modules(logfile)
					results = step.run(sample, logfile)
					if results is not None and isinstance(results, dict):
						for label, path in results.iteritems():
							sample.add_file(step.name, label, path)
				except Exception as e:
					if step.is_critical():
						status_counts['critical'] += 1
						raise #this step was critical and failed, rethrow to the outer try and halt the pipeline!
					else:
						status_counts['warn'] += 1
						logfile.write('Encountered error during processing:\n')
						logfile.write('%s\n' % e)
						logfile.flush()
				finally:
					logfile.write("--------------------------------------------\n\n")
					logfile.flush()

		except Exception as e:
			tasks_statuses[sample.name] = tasks_statuses[sample.name].update(None, 'Error!').finish()
			logfile.write('Encountered error during processing:\n')
			logfile.write('%s\n' % e)
			logfile.write(traceback.format_exc())
		else:
			#let the user know our progress
			tasks_statuses[sample.name] = tasks_statuses[sample.name].update(1, 'Done!').finish()
			logfile.write('Completed processing of sample "%s"!\n' % (sample.name,))
			logfile.write('-> See output: "%s"\n' % (sample.dest,))
		finally:
			logfile.write("-> Wall clock: %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"),))
			logfile.write("--------------------------------------------\n")
			logfile.write("Total Pipeline Steps: %d\n" % (status_counts['total'],))
			logfile.write("-> # Steps Run:  %d\n" % (status_counts['attempted'],))
			logfile.write("-> # Steps Warn: %d (non-critical failure)\n" % (status_counts['warn'],))
			logfile.write("-> # Steps Fail: %d (critical failure)\n" % (status_counts['critical'],))
			logfile.write("--------------------------------------------\n\n")
			logfile.flush()
			return sample
#end __execute_pipeline_on_sample()
