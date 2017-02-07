#!/usr/bin/python
# -*- coding: UTF-8 -*-
import os
import sys
import Common
import time
import platform
import traceback
import re
#from multiprocessing import Pool, Manager 
#from multiprocessing import Manager
import dill	#use dill for pickling, actually supports serializing useful things! (i.e. lambdas, objects)
import multiprocess as mp	#use this ls alternative multiprocessing from pathos, used in combination with dill
from ProcessHelpers import MultiStatusProgressBar, MultiStatusProgressItem
import subprocess


_global_manager = mp.Manager()


class GenomeInfo:
	def __init__(self, name, gsize, chrsize=None):
		'''genome name, typically UCSC naming convention '''
		self.name = name
		'''Effective genome size'''
		self.gsize = gsize
		'''location of chromosome size info file'''
		self.chrsize = chrsize
		self.indexes = {}
		self.wg_fasta = None
		self.chr_fasta = {}
	#end __init__()
	
	def try_discover(self, basepath):
		import glob
		'''given a base path, try to discover indexes and other reference genome data according to the illumina golden path layout'''
		#find indexes
		idx_types = {
			'BowtieIndex': 	'*.1.ebwt',
			'Bowtie2Index':	'*.1.bt2', 
			'BWAIndex': 	'*.fa.bwt'
		}
		for idx in idx_types.keys():
			idx_path = os.path.join(basepath, 'Sequence', idx)
			if os.path.exists(idx_path):
				matches = glob.glob(os.path.join(idx_path, idx_types[idx]))
				if len(matches) > 0:
					self.add_index(idx, os.path.join(idx_path, Common.basename_noext(matches[0], True)))
		
		#chromosome fasta files
		for f in glob.glob(os.path.join(basepath, 'Sequence', 'Chromosomes', '*.fa')):
			self.chr_fasta[Common.basename_noext(f)] = f
		
		#whole genome fasta
		wg_fasta_results = glob.glob(os.path.join(basepath, 'Sequence', 'WholeGenomeFasta', '*.fa'))
		if len(wg_fasta_results) > 0:
			self.wg_fasta = wg_fasta_results[0]
		
		#chromsizes
		if self.chrsize is not None:
			chrsize_results = glob.glob(os.path.join(basepath, 'Sequence', 'WholeGenomeFasta', 'chrom.sizes'))
			if len(chrsize_results) > 0:
				self.chrsize = chrsize_results[0]
	#end try_discover()
	
	def add_index(self, name, value):
		self.indexes[name] = value
	#end add_index
	
	def has_index(self, name):
		return name in self.indexes
	#end has_index()
	
	def get_index(self, name):
		if self.has_index(name):
			return self.indexes[name]
		return None
	#end get_index()
	
	__known_references = None
	__probe_paths = []
	@staticmethod
	def register_path(path):
		GenomeInfo.__probe_paths.append(path)
	#end register_path()
	
	@staticmethod
	def get_reference_genomes():
		if GenomeInfo.__known_references is None:
			GenomeInfo.__known_references = {}
			
			conanocal = []
			conanocal.append(GenomeInfo('hg18', '2.77e9'))
			conanocal.append(GenomeInfo('hg19', '2.79e9'))
			conanocal.append(GenomeInfo('mm9',  '1.91e9'))
			conanocal.append(GenomeInfo('mm10', '2.15e9'))
			
			for c in conanocal:
				for p in GenomeInfo.__probe_paths:
					#print "trying path %s for genome %s" % (os.path.join(p, c.name), c.name)
					if os.path.exists(os.path.join(p, c.name)):
						#print "Found!"
						c.try_discover(os.path.join(p, c.name))
						GenomeInfo.__known_references[c.name] = c
						break
		return GenomeInfo.__known_references
	#end get_reference_genomes()
#end class GenomeInfo
GenomeInfo.register_path('/mnt/ref/reference/Homo_sapiens/UCSC')
GenomeInfo.register_path('/mnt/ref/reference/Mus_musculus/UCSC')
GenomeInfo.register_path('/mnt/ref/reference/Rattus_norvegicus/UCSC')
GenomeInfo.register_path('/mnt/ref/reference/Saccharomyces_cerevisiae/UCSC')
GenomeInfo.register_path('/mnt/ref/reference/PhiX/UCSC')

GenomeInfo.register_path('/home/thackray/reference/Homo_sapiens/UCSC')
GenomeInfo.register_path('/home/thackray/reference/Mus_musculus/UCSC')
GenomeInfo.register_path('/home/thackray/reference/Rattus_norvegicus/UCSC')
GenomeInfo.register_path('/home/thackray/reference/Saccharomyces_cerevisiae/UCSC')
GenomeInfo.register_path('/home/thackray/reference/PhiX/UCSC')



class FileInfo:

	def __init__(self, filepath):
		self._fullpath = filepath
		self._companions = _global_manager.list()
	
	@property
	def fullpath(self):
		return self._fullpath
	
	@property
	def dirname(self):
		return os.path.dirname(self._fullpath)
	
	@property
	def basename(self):
		return os.path.basename(self._fullpath)
		
	@property
	def filename(self):
		return os.path.splitext(self.basename)[0]
		
	@property
	def ext(self):
		return os.path.splitext(self.basename)[1]
		
	def has_companions(self):
		return len(self._companions) > 0
		
	def add_companion(self, filepath):
		self._companions.append(FileInfo(filepath))
#end class FileInfo
		

class PipelineSample:

	def __init__(self, name, genome, dest, format):
		self.files = _global_manager.dict()
		self.attr = _global_manager.dict()
		self.name = name
		if isinstance(genome, GenomeInfo):
			self.genome = genome
		else:
			self.genome = GenomeInfo.get_reference_genomes()[genome]
		self.dest = os.path.abspath(dest)
		self.format = format.lower()
	#end __init__()
	
	def __getstate__(self):
		state = dict(self.__dict__)
		state['files'] = dict(self.files)
		state['attr'] = dict(self.attr)
		return state
	#end __getstate__()
	
	def __setstate__(self, state):
		self.files = _global_manager.dict(state['files'])
		self.attr = _global_manager.dict(state['attr'])
		self.__dict__.update(state)
	#end __setstate__()
	
	def add_file(self, group, label, path):
		if group not in self.files:
			self.files[group] = {}
		#necessary to swap the whole sub-dict when multiprocessing!
		tmp = self.files[group]
		tmp[label] = path
		self.files[group] = tmp
	#end add_file()
	
	def has_file_group(self, group):
		return group in self.files
	#end has_file_group()
	
	def has_file(self, group, label):
		if group in self.files:
			if label in self.files[group]:
				return True
		return False
	#end has_file()
	
	def get_file(self, group, label):
		if self.has_file(group, label):
			return self.files[group][label] 
		return None
	#end has_file()
	
	def get_file_group(self, group):
		if group in self.files:
			return self.files[group]
		return None
	#end get_file_group()
	
	def add_attribute(self, name, value):
		self.attr[name] = value
	#end add_attribute()
	
	def has_attribute(self, name):
		return name in self.attr
	#end has_attribute()
	
	def get_attribute(self, name):
		if self.has_attribute(name):
			return self.attr[name]
		return None
	#end get_attribute()
	
	def remove_attribute(self, name):
		if self.has_attribute(name):
			del self.attr[name]
		return None
	#end remove_attribute()
#end PipelineSample

class ModuleParameter:
	def __init__(self, name, type, default, value=None, desc="", nullable=False):
		self.name = name
		self.type = type
		self.description = desc
		self.nullable = nullable
		self.default = self._coerce_value(default, self.type)
		if value is None and not self.nullable:
			self.value = self.default
		else:
			self.value = value
	#end __init__()
	
	def is_defualt(self):
		return default == self._coerce_value(self.value)
	#end is_defualt()
	
	def reset(self):
		self.value = self.default
	#end is_defualt()
	
	def get_value(self):
		return self._coerce_value(self.value, self.type)
	#end get_value()
	
	def get_value_as_type(self, type):
		return self._coerce_value(self.get_value(), type)
	#end get_value()

	def _coerce_value(self, value, type):
		if self.nullable and value is None:
			return None
		elif type == bool:
			from distutils.util import strtobool
			return strtobool(str(value))
		elif type == int:
			return int(value)
		elif type == float:
			return float(value)
		elif type == long:
			return long(value)
		elif type == complex:
			return complex(value)
		elif type == str:
			return str(value)
		elif type == list:
			return list(value)
		elif type == tuple:
			return tuple(value)
		elif type == set:
			return set(value)
		else:
			raise ValueError('Unable to convert %s to %s for parameter %s!' % (str(value), str(type), str(self.name),))
	#end _coerce_value()
	
	def __str__(self):
		return "(%s)%s = %s [%s] %s" % (str(self.type), str(self.name), str(self.get_value()), str(self.default), str(self.description))
	#end __str__()
	
	def __repr__(self):
		return "ModuleParameter(%s, %s, %s, %s, %s)" % (str(self.name), str(self.type), str(self.default), str(self.get_value()), str(self.description))
	#end __repr__()
#end class ModuleParameter


class PipelineModule:

	def __init__(self, name, short_description):
		self.name = name
		self.description = short_description
		self.processors = 1
		self.parameters = {}
		self.resolvers = {}
		self._critical = False
	#end __init__()
	
	def set_critical(self, is_critical):
		self._critical = is_critical
	#end set_critical()
	
	def is_critical(self):
		return self._critical
	#end is_critical()
	
	def add_parameter(self, parameter):
		self.parameters[parameter.name] = parameter
	#end add_parameter()
	
	def set_parameter(self, parameter_name, value):
		if self.has_parameter(parameter_name):
			self.parameters[parameter_name].value = value
		else:
			raise ValueError('Parameter "%s" does not exist!' % (parameter_name,))
	#end add_parameter()
	
	def has_parameter(self, parameter_name):
		return parameter_name in self.parameters
	#end has_parameter()
	
	def get_parameter(self, parameter_name):
		if self.has_parameter(parameter_name):
			return self.parameters[parameter_name]
		else:
			raise ValueError('Parameter "%s" does not exist!' % (parameter_name,))
	#end get_parameter()
	
	def get_parameter_value(self, parameter_name):
		if self.has_parameter(parameter_name):
			return self.parameters[parameter_name].get_value()
		else:
			raise ValueError('Parameter "%s" does not exist!' % (parameter_name,))
	#end get_parameter()
	
	def get_parameter_value_as_string(self, parameter_name):
		if self.has_parameter(parameter_name):
			return self.parameters[parameter_name].get_value_as_type(str)
		else:
			raise ValueError('Parameter "%s" does not exist!' % (parameter_name,))
	#end get_parameter()
	
	def _name_resolver(self, name):
		self.resolvers[name] = lambda sample: None
	#end _name_resolver()
	
	def set_resolver(self, name, function):
		self.resolvers[name] = function
	#end set_resolver()
	
	def resolve_input(self, name, sample):
		return self.resolvers[name](sample)
	#end resolve_input()
	
	def get_resolver_names(self):
		return self.resolvers.keys()
	#end get_resolver_names()
	
	def set_available_cpus(self, cpus):
		self.processors = cpus
	#end set_available_processors()
	
	def show_version(self, handle=None, fancy=True):
		pass
	#end show_version()
	
	def load_modules(self):
		""" Loads any system modules required for this module to function. """
		pass
	#end load_modules()
	
	def run(self, sample, logfile):
		""" Runs the pipeline module on the provided PipelineSample and sending log information to logfile handle """
		pass
	#end run()
	
	def supported_types(self):
		""" Returns an iterable of the file types this module supports """
		pass
	#end supported_types()
	
	def is_compatible_with_sample(self, sample):
		#if self.supported_types() is not None:
		#	return sample.format in self.supported_types()
		return True
	#end is_compatible_with_sample()
	
	def _run_subprocess(self, cmd, **kwargs):
		proc = subprocess.Popen(cmd, **kwargs)
		proc.communicate()
		
		if not proc.returncode == 0:
			raise subprocess.CalledProcessError(proc.returncode, str(cmd))
	#end _run_subprocess()
	
	def documentation(self):
		hash_length = 40
		buffer = "%s\n%s\n%s\n%s\n" % (self.name, '-'*hash_length, self.description, '-'*hash_length)
		buffer += "Parameters:"
		if len(self.parameters) < 1:
			buffer += "\n\tNo Parameters Declared\n"
		else:
			buffer += "\n"
			for param in self.parameters:
				buffer += "\t%s: %s\n" % (param, str(self.parameters[param]))
		#buffer += '-'*hash_length
		buffer += '\n'
		
		buffer += "Resolvers:"
		if len(self.resolvers) < 1:
			buffer += "\n\tNo Resolvers Declared\n"
		else:
			buffer += "\n"
			for resolver in self.resolvers:
				buffer += "\t%s: %s\n" % (resolver, str(self.resolvers[resolver]))
		buffer += '-'*hash_length
		buffer += '\n'
		return buffer
	#end documentation()
#end class PipelineModule



class AnalysisPipeline:

	def __init__(self, name="Anonymous Pipeline"):
		self.name = name
		self.pipeline = []
		self.tasks_statuses = None
	#end __inti__()
	
	def safe_name(self):
		return re.sub('[^\w\-_\.]', '_', self.name)
	#end safe_name()
	
	def append_module(self, module, critical=False):
		module.set_critical(critical)
		self.pipeline.append(module)
	#end append_module()
	
	def insert_module(self, position, module, critical=False):
		module.set_critical(critical)
		self.pipeline.insert(position, module)
	#end append_module()
	
	def documentation(self):
		hash_length = 40
		buffer = "%s\nBegin Pipeline: %s\n%s\n" % ('='*hash_length, self.name, '='*hash_length)
		
		for m in self.pipeline:
			buffer += "\n%s\n\n" % (u'↓'*hash_length,)
			buffer += m.documentation()
		
		buffer += "\n%s\n\n%s\nEnd Pipeline: %s\n%s\n" % (u'↓'*hash_length,'='*hash_length, self.name, '='*hash_length)
		return buffer
	#end explain()
#end class AnalysisPipeline

class PipelineRunner:
	def __init__(self, pipeline):
		self.pipeline = pipeline
		self.tasks_statuses = None

	def run(self, samples):
		pass
#end class PipelineRunner

class SerialPipelineRunner(PipelineRunner):
	def __init__(self, pipeline):
		PipelineRunner.__init__(self, pipeline)

	def run(self, samples):
		for sample in samples:
			_execute_pipeline_on_sample(self.pipeline, sample, {sample.name: MultiStatusProgressItem(sample.name, 'Queued...')})
#end class SerialPipelineRunner

class ParallelPipelineRunner(PipelineRunner):
	def __init__(self, pipeline, threads):
		PipelineRunner.__init__(self, pipeline)
		self.threads = threads

	def run(self, samples):
		sample_count = len(samples)
		if sample_count < 1:
			return #No samples to run!
		
		self.tasks_statuses = _global_manager.dict()
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
#end class ParallelPipelineRunner

class SlurmPipelineRunner(PipelineRunner):
	def __init__(self, pipeline):
		PipelineRunner.__init__(self, pipeline)

	def run(self, samples, nodes, threads_per_node):
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
			self.tasks_statuses = _global_manager.dict()
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
					'--partition', 'main',
					'-n', '1',
					'--cpus-per-task', str(threads_per_node),
					#'--export', 'ALL',
					'python', os.path.abspath(__file__),
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
#end class SerialPipelineRunner


def _execute_pipeline_on_sample(pipeline, sample, tasks_statuses):
	tasks_statuses[sample.name] = tasks_statuses[sample.name].start()
	#ensure the output directory exists
	Common.ensure_dir(sample.dest)

	with open(os.path.join(sample.dest, sample.name+'.log'), 'a', buffering=0) as logfile:
		sys.stdout = sys.stderr = logfile
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
					step.load_modules()
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


if __name__ == "__main__":
	import argparse
	import threading
	
	parser = argparse.ArgumentParser()
	parser.add_argument('pipeline', help="pipeline pickle file")
	parser.add_argument('sample', help="sample pickle file")
	parser.add_argument('status', help="sample status pickle file")
	args = parser.parse_args()
	
	if not os.path.exists(args.pipeline):
		raise ValueError("Path %s not found!" % (args.pipeline,))
	if not os.path.exists(args.sample):
		raise ValueError("Path %s not found!" % (args.sample,))
	if not os.path.exists(args.status):
		raise ValueError("Path %s not found!" % (args.status,))
	
	with open(args.pipeline, 'rb') as pf:
		pipeline_pickle = dill.load(pf)
	with open(args.sample, 'rb') as sf:
		sample_pickle = dill.load(sf)
	with open(args.status, 'rb') as ssf:
		status_pickle = dill.load(ssf)
	
	tasks_statuses = dict()
	tasks_statuses[sample_pickle.name] = status_pickle
	
	thread = threading.Thread(target=_execute_pipeline_on_sample, args=(pipeline_pickle, sample_pickle, tasks_statuses))
	thread.start()
	
	while thread.is_alive():
		time.sleep(0.5)
		with open(args.status, 'wb', 0) as f:
			dill.dump(tasks_statuses[sample_pickle.name], f)
	thread.join()
	

