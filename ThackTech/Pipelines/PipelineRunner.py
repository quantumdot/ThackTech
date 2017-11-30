import os
import sys
import traceback
import platform
import time
from ThackTech import filetools, conf
from ThackTech.Pipelines import AnalysisPipeline, PipelineSample, ModuleRunContext, FileInfo, FileContext, CPU_COUNT


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
	
	def monitor(self):
		pass
#end class PipelineRunner



def add_runner_args(argparser):
	"""Adds an argument group to the passed argparse object with pipeline runner options
	
	Parameters:
		argparser:	An argparse.ArgumentParser object to add runner options to
		
	Returns:
		argument group for performance options
	"""
	default_args = {
		'general': {
			'runner': 'serial',
			'shm_dir': '/run/shm'
		},
		'slurm_runner': {
			'partition': "main",
			'nodes': 1,
			'threads': CPU_COUNT,
			'time_limit': "1:00:00",
			'mem': '4G',
			'cmd': 'sbatch'
		},
		'parallel_runner': {
			'threads': CPU_COUNT
		},
		'serial_runner': {
			'threads': CPU_COUNT
		}
	}
	
	
	cfg = conf.get_config('pipelines')
	for section in list(default_args.keys()):
		for option in list(default_args[section].keys()):
			try:
				default_args[section][option] = cfg.get(section, option)
			except:
				pass
	
	runner_type = default_args['general']['runner']+'_runner'
	default_threads = default_args[runner_type]['threads']
		
	performance_group = argparser.add_argument_group('Performance')
	performance_group.add_argument('-p', '--threads', type=int, default=default_threads, help="Number of processors to use for processing.")
	performance_group.add_argument('--shm', action='store_true', help="Use ramfs for file IO.")
	performance_group.add_argument('--shm-path', action='store', default=default_args['general']['shm_dir'], help='When --shm is passed, the path to use for ram disk storage. Individual samples will have dedicated subfolders on this path. Please ensure this path has appropriate permissions.')
	performance_group.add_argument('--runner', action='store', default=default_args['general']['runner'], choices=['slurm', 'parallel', 'serial'], help="Which pipeline runner to use.")
	
	performance_group.add_argument('--slurm-partition', action='store', default=default_args['slurm_runner']['partition'], help="For slurm runner, the partition to run jobs on.")
	performance_group.add_argument('--slurm-time', action='store', default=default_args['slurm_runner']['time_limit'], help="For slurm runner, time limit for jobs.")
	performance_group.add_argument('--slurm-mem', action='store', default=default_args['slurm_runner']['mem'], help="For slurm runner, memory for jobs.")
	performance_group.add_argument('--slurm-cmd', action='store', default=default_args['slurm_runner']['cmd'], help="For slurm runner, style of job invocation. One of [srun, sbatch].")
	
	
	performance_group.add_argument('--module-config', action='append', default=[], help="Specify additional pipeline module config files.")
	performance_group.add_argument('--pipeline-config', action='append', default=[], help="Specify additional pipeline config files.")
	performance_group.add_argument('--genome-config', action='append', default=[], help="Specify additional genome config files.")
	return performance_group
#end get_runner_args()


def get_configured_runner(args, pipeline, **kwargs):
	
	if args.runner == 'slurm':
		from ThackTech.Pipelines import SlurmPipelineRunner
		runner = SlurmPipelineRunner(pipeline, partition=args.slurm_partition, nodes=1, threads=args.threads, time_limit=args.slurm_time, mem=args.slurm_mem, slurm_cmd=args.slurm_cmd, **kwargs)
	else:
		do_slurm_safety_check()
		
		if args.runner == 'parallel':
			from ThackTech.Pipelines import ParallelPipelineRunner
			runner = ParallelPipelineRunner(pipeline, args.threads)
		else: #serial runner
			from ThackTech.Pipelines import SerialPipelineRunner
			runner = SerialPipelineRunner(pipeline)
	
	return runner
#end get_configured_runner()

def do_slurm_safety_check():
	#call some slurm command to check if slurm present
	
	#if slurm present:
	#	alert user
	#	ask for confirmation that you want to use non-slurm runner on a slurm system
	#	if user confirms usage:
	#		continue
	#	otherwise:
	#		die
	pass
#end do_slurm_safety_check()



def _execute_pipeline_on_sample(pipeline, sample, tasks_statuses):
	"""Exexutes a pipeline on a given sample
	
	This method is responsible for directly running a pipeline on a given sample, although
	it is not intended to use this method directly. Instead, one should use a derivitive 
	of the PipelineRunner class (e.x. SerialPipelineRunner, ParallelPipelineRunner, or 
	SlurmPipelineRunner).
	
	Parameters
		pipeline:	(AnalysisPipeline) Pipeline to run
		sample:		(PipelineSample) sample to run pipeline on
		task_statuses: (complicated....) Status information
	"""
	#some assertions to make sure we don't go too crazy!
	assert isinstance(pipeline, AnalysisPipeline), "pipeline parameter must be of type AnalysisPipeline!"
	assert isinstance(sample, PipelineSample), "sample parameter must be of type Pipelinesample!"
	
	#ensure the output directory exists
	filetools.ensure_dir(sample.dest)

	with open(os.path.join(sample.dest, sample.name+'.log'), 'a', buffering=0) as logfile:
		sys.stdout = sys.stderr = logfile
		try:
			tasks_statuses[sample.name] = tasks_statuses[sample.name].start()
			output_manifest_location = os.path.join(sample.dest, sample.name+'_output_manifest.tsv')
			logfile.write(pipeline.documentation())
			
			pipeline_size = len(pipeline)
			pipeline_steps = pipeline.itersteps()
			status_counts = {
				'total': 	pipeline_size,
				'skipped':  0,
				'attempted':0,
				'warn':		0,
				'critical':	0
			}
		
			tasks_statuses[sample.name] = tasks_statuses[sample.name].update(0, 'Preparing...')
			logfile.write('Processing sample "{}"....\n'.format(sample.name))
			logfile.write('-> Pipeline: {} ({} steps)\n'.format(pipeline.name, pipeline_size))
			
			logfile.write("{}".format(pipeline_steps))
				
			if pipeline.offset is not None:
				status_counts['skipped'] = pipeline_steps[0].index
				logfile.write('-> Resuming from step: {} (checkpoint "{}")\n'.format(pipeline_steps[0].index+1, pipeline.offset))
				sample.read_file_manifest(output_manifest_location)
			
			logfile.write('-> Running on: {}\n'.format(' '.join(platform.uname())))
			logfile.write("-> Wall clock: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
			logfile.write("--------------------------------------------\n\n")
			
				
			
			for step in pipeline_steps:
				status_counts['attempted'] += 1
				logfile.write('Running pipeline step #{}: {}\n'.format(step.index+1, step.module.name))
				logfile.write("-> Wall clock: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
				logfile.flush()
				try:
					#if not step.is_compatible_with_sample(sample):
					#	raise ValueError('%s does not support the "%s" file format!\nThe following types ARE supported: %s\nSkipping processing......\n' % (step.name, sample.format, ', '.join(step.supported_types())))
				
					tasks_statuses[sample.name] = tasks_statuses[sample.name].update(float(step.index)/float(pipeline_size), step.module.description)
					#step.load_modules(logfile)
					cxt = ModuleRunContext(pipeline.name, step.index, step.module.name, logfile, sample)
					results = step.module.run(cxt)
					
					if results is not None:
						
						if isinstance(results, dict):
							#support the old-style of returning output
							#this may be removed in the future!
							for label, path in results.items():
								sample.add_file(FileInfo(path, FileContext.from_module_context(cxt, label)))
						else:
							for f in results:
								sample.add_file(f)
						
						#done running this step, update the output manifest
						sample.write_file_manifest(output_manifest_location)
			
				except Exception as e:
					if step.module.is_critical:
						status_counts['critical'] += 1
						raise #this step was critical and failed, rethrow to the outer try and halt the pipeline!
					else:
						status_counts['warn'] += 1
						logfile.write('Encountered error during processing:\n')
						logfile.write('{}\n'.format(e))
						logfile.write(traceback.format_exc())
						logfile.write('\n\nModule is not marked critical, continuing with pipeline....\n')
						logfile.flush()
				finally:
					logfile.write("--------------------------------------------\n\n")
					logfile.flush()
					
		except Exception as e:
			tasks_statuses[sample.name] = tasks_statuses[sample.name].update(None, 'Error!').finish()
			logfile.write('Encountered error during processing:\n')
			logfile.write('{}\n'.format(e))
			logfile.write(traceback.format_exc())
			logfile.write('\n')
		else:
			#let the user know our progress
			tasks_statuses[sample.name] = tasks_statuses[sample.name].update(1, 'Done!').finish()
			logfile.write('Completed processing of sample "{}"!\n'.format(sample.name))
			logfile.write('-> See output: "{}"\n'.format(sample.dest))
		finally:
			logfile.write("-> Wall clock: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
			logfile.write("--------------------------------------------\n")
			logfile.write("Total Pipeline Steps: {}\n".format(status_counts['total']))
			logfile.write("-> # Steps Skipped:  {}\n".format(status_counts['skipped']))
			logfile.write("-> # Steps Run:  {}\n".format(status_counts['attempted']))
			logfile.write("-> # Steps Warn: {} (non-critical failure)\n".format(status_counts['warn']))
			logfile.write("-> # Steps Fail: {} (critical failure)\n".format(status_counts['critical']))
			logfile.write("--------------------------------------------\n\n")
			logfile.flush()
			return sample
#end __execute_pipeline_on_sample()
