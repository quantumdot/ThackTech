import os
import sys
import traceback
import platform
import time
import dill	#use dill for pickling, actually supports serializing useful things! (i.e. lambdas, objects)
import multiprocess as mp	#use this ls alternative multiprocessing from pathos, used in combination with dill
from ThackTech import Common

# try:
	# #try to add module system to python path
	# sys.path.insert(0, os.path.join(os.environ['MODULESHOME'], "init"))
# except:
	# pass

class PipelineRunner:
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
	Common.ensure_dir(sample.dest)

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


# if __name__ == "__main__":
	# import argparse
	# import threading
	# import dill
	
	# parser = argparse.ArgumentParser()
	# parser.add_argument('pipeline', help="pipeline pickle file")
	# parser.add_argument('sample', help="sample pickle file")
	# parser.add_argument('status', help="sample status pickle file")
	# args = parser.parse_args()
	
	# if not os.path.exists(args.pipeline):
		# raise ValueError("Path %s not found!" % (args.pipeline,))
	# if not os.path.exists(args.sample):
		# raise ValueError("Path %s not found!" % (args.sample,))
	# if not os.path.exists(args.status):
		# raise ValueError("Path %s not found!" % (args.status,))
	
	# with open(args.pipeline, 'rb') as pf:
		# pipeline_pickle = dill.load(pf)
	# with open(args.sample, 'rb') as sf:
		# sample_pickle = dill.load(sf)
	# with open(args.status, 'rb') as ssf:
		# status_pickle = dill.load(ssf)
	
	# tasks_statuses = dict()
	# tasks_statuses[sample_pickle.name] = status_pickle
	
	# thread = threading.Thread(target=_execute_pipeline_on_sample, args=(pipeline_pickle, sample_pickle, tasks_statuses))
	# thread.start()
	
	# while thread.is_alive():
		# time.sleep(0.5)
		# with open(args.status, 'wb', 0) as f:
			# dill.dump(tasks_statuses[sample_pickle.name], f)
	# thread.join()
	

