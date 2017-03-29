import subprocess
from tabulate import tabulate
from ThackTech.Pipelines import ModuleParameter
from ThackTech import conf
import textwrap
from abc import ABCMeta, abstractmethod
from _pyio import __metaclass__


class PipelineModule(object):
	"""Represents a self-contained work unit within a pipeline
	
	
	
	"""
	__metaclass__ = ABCMeta #enable abstract method decarator
	
	def __init__(self, name, short_description, critical=False, processors=1):
		self.name = name
		self.description = short_description
		self.parameters = {}
		self.resolvers = {}
		self.processors = processors
		self._critical = critical
	#end __init__()
	
	
	
	def tool_versions(self):
		"""Return a dict of tool versions used by this module
		
		dictionary is of the form <ToolName> -> <VersionString>
		"""
		return {}
	#end show_version()
	
	def load_modules(self, logfile):
		""" Loads any system modules required for this module to function. """
		pass
	#end load_modules()
	
	@abstractmethod
	def run(self, context):
		""" Runs the pipeline module given the supplied context
		
		!!!!! This method is the most important method to override in concrete implementations !!!!!
		
		Context provide particular run specific information, including:
			-sample (PipelineSample) and sending log information to logfile handle 
			-logfile (file-like)
			-misc. metadata about step number, pipeline name, etc.
			
		Parameters:
			context: (ModuleRunContext) Context to run this module with.
			
		Returns:
			iterable of FileInfo objects representing output files
		"""
		pass
	#end run()
	
	
	
	
	
	
	
	
	
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
	
	def set_parameters_from_config(self):
		"""Reads from the pipeline_modules set of configuration files and sets this modules parameters accordingly
		
		Configuration data is interpreted from config files in the following way:
			files named "pipeline_modules[.default].ini
			config section name corresponds to module name
		""" 
		cfg = conf.get_config('pipeline_modules')
		if cfg.has_section(self.name):
			for param in self.parameters.itervalues():
				if cfg.has_option(self.name, param.name):
					self.set_parameter(param.name, cfg.get(self.name, param.name))
	#end set_parameters_from_config()
	
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
	
	def resolve_input(self, name, cxt):
		return self.resolvers[name](cxt)
	#end resolve_input()
	
	def get_resolver_names(self):
		"""Gets a list of declared resolver names
		"""
		return self.resolvers.keys()
	#end get_resolver_names()
	
	
	
	
	
	def set_available_cpus(self, cpus):
		self.processors = cpus
	#end set_available_processors()
	
	
	
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
		buff = ""
		hash_length = 40
		
		buff += "%s\n" % (self.name,)
		buff += "%s\n" % ('-'*hash_length,)
		buff += "CRITICAL:    %s\n" % (self.is_critical,)
		buff += "PROCESSORS:  %d\n" % (self.processors) 
		buff += "%s\n"  % (textwrap.fill("DESCRIPTION: "+self.description, hash_length),)
		buff += "%s\n" % ('-'*hash_length,)
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
		buff += '\n'
		
		buff += "TOOL VERSIONS:"
		tools = self.tool_versions()
		if len(self.resolvers) < 1:
			buff += "\n\tNo Tools Declared\n"
		else:
			buff += "\n"
			for tool, version in tools:
				buff += "\t%s: %s\n" % (tool, version)
		buff += '-' * hash_length
		buff += '\n'
		return buff
	#end documentation()
#end class PipelineModule