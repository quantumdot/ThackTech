import os
import subprocess
from tabulate import tabulate
from ThackTech.Pipelines import ModuleParameter


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
	
	def load_modules(self, logfile):
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
		#print os.environ
		proc = subprocess.Popen(cmd, **kwargs)
		out, err = proc.communicate()
		
		if not proc.returncode == 0:
			raise subprocess.CalledProcessError(proc.returncode, str(cmd), str(out)+"\n"+str(err))
		return (out, err)
	#end _run_subprocess()
	
	def documentation(self):
		hash_length = 40
		buffer = "%s\n%s\n%s\n%s\n" % (self.name, '-'*hash_length, self.description, '-'*hash_length)
		buffer += "PARAMETERS:"
		if len(self.parameters) < 1:
			buffer += "\n\tNo Parameters Declared\n"
		else:
			buffer += "\n"
			param_headers = ['Name', 'Type', 'Value', 'Default', 'Description']
			param_table = []
			for param in self.parameters.itervalues():
				param_table.append([str(param.name), param.get_type_as_string(), str(param.get_value()), str(param.default), str(param.description)])
				#buffer += "\t%s: %s\n" % (param, str(self.parameters[param]))
			buffer += tabulate(param_table, headers=param_headers, tablefmt="simple")
			buffer += '\n'
		#buffer += '-'*hash_length
		buffer += '\n'
		
		buffer += "RESOLVERS:"
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