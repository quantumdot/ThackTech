import subprocess
from tabulate import tabulate
from ThackTech.Pipelines import ModuleParameter


class PipelineModule(object):
	"""Represents a self-contained work unit within a pipeline
	
	
	
	"""
	def __init__(self, name, short_description):
		self.name = name
		self.description = short_description
		self.processors = 1
		self.parameters = {}
		self.resolvers = {}
		self._critical = False
	#end __init__()
	
	def set_critical(self, is_critical):
		"""Sets if failure of this module is a critical failure.
		If this module fails and is marked as critical, then the remainder of
		the pipeline will be aborted. Useful if downstream steps rely on the 
		output of this module.
		"""
		self._critical = is_critical
	#end set_critical()
	
	@property
	def is_critical(self):
		"""Returns bool where True indicates that this module is critical
		"""
		return self._critical
	#end is_critical()
	
	def add_parameter(self, parameter):
		"""Adds a parameter to this module
		Typically this is used only derivitives of PipelineModule.
		To set the value of a parameter, use the set_parameter() method
		"""
		if not isinstance(parameter, ModuleParameter):
			raise ValueError('Expecting parameter to be of type ModuleParmater, %s given!' % (type(parameter).__name__,))
		self.parameters[parameter.name] = parameter
	#end add_parameter()
	
	def set_parameter(self, parameter_name, value):
		if self.has_parameter(parameter_name):
			self.parameters[parameter_name].value = value
		else:
			raise ValueError('Parameter "%s" does not exist!' % (parameter_name,))
	#end set_parameter()
	
	def has_parameter(self, parameter_name):
		"""Check if a parameter by the name of `parameter_name` has been declared by this module
		"""
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
			return self.parameters[parameter_name].value
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
		"""Names a resolver that should be supplied by the consumer
		
		When a resolver is initially named, the defualt implementation
		assigned to the resolver is a lambda that returns None.
		"""
		self.resolvers[name] = lambda sample: None
	#end _name_resolver()
	
	def set_resolver(self, name, _callable):
		self.resolvers[name] = _callable
	#end set_resolver()
	
	def resolve_input(self, name, sample):
		return self.resolvers[name](sample)
	#end resolve_input()
	
	def get_resolver_names(self):
		"""Gets a list of declared resolver names
		"""
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
		"""Deprecated
		"""
		#if self.supported_types() is not None:
		#	return sample.format in self.supported_types()
		return True
	#end is_compatible_with_sample()
	
	def _run_subprocess(self, cmd, **kwargs):
		"""Wrapper around the subprocess.Popen call and provides some extra convience
		
		Use is basically the same as the subprocess.Popen call, except that this method
		checks the return code of the completed process and raises a CalledProcessError
		if the return code <> 0. 
		
		This method returns a tuple of (stdout, stderr) from the process, but this requires
		setting stdout and stderr to subprocess.PIPE in kwargs.
		
		This method will block until the process has completed!
		"""
		#print os.environ
		proc = subprocess.Popen(cmd, **kwargs)
		out, err = proc.communicate()
		
		if not proc.returncode == 0:
			raise subprocess.CalledProcessError(proc.returncode, str(cmd), str(out)+"\n"+str(err))
		return (out, err)
	#end _run_subprocess()
	
	def documentation(self):
		"""Return a string that documents this module and its parameters and resolvers
		"""
		hash_length = 40
		buff = "%s\n%s\n%s\n%s\n" % (self.name, '-'*hash_length, self.description, '-'*hash_length)
		buff += "PARAMETERS:"
		if len(self.parameters) < 1:
			buff += "\n\tNo Parameters Declared\n"
		else:
			buff += "\n"
			param_headers = ['Name', 'Type', 'Value', 'Default', 'Nullable', 'Choices', 'Description', ]
			param_table = []
			for param in self.parameters.itervalues():
				param_table.append([str(param.name), param.type_name, str(param.value), str(param.default), str(param.nullable), str(param.choices), str(param.description)])
				#buffer += "\t%s: %s\n" % (param, str(self.parameters[param]))
			buff += tabulate(param_table, headers=param_headers, tablefmt="simple")
			buff += '\n'
		#buffer += '-'*hash_length
		buff += '\n'
		
		buff += "RESOLVERS:"
		if len(self.resolvers) < 1:
			buff += "\n\tNo Resolvers Declared\n"
		else:
			buff += "\n"
			for resolver in self.resolvers:
				buff += "\t%s: %s\n" % (resolver, str(self.resolvers[resolver]))
		buff += '-' * hash_length
		buff += '\n'
		return buff
	#end documentation()
#end class PipelineModule