from multiprocess.sharedctypes import _new_value
class ModuleParameter(object):
	"""Represents a parameter for a PipelineModule
	
	Provides type-safety, default, nullability, and parameter documentation
	"""
	def __init__(self, name, value_type, default, value=None, desc="", nullable=False, choices=None):
		self.name = name
		self.type = value_type
		self.description = desc
		self.nullable = nullable
		self.choices = choices
		self.default = self._coerce_value(default, self.type)
		if value is None and not self.nullable:
			self.value = self.default
		else:
			self.value = value
	#end __init__()
	
	@property
	def value(self):
		return self.__value
	
	@value.setter
	def value(self, value):
		self.__value = self._coerce_value(value, self.type)

	@value.deleter
	def value(self):
		self.value = self.default
	
	@property
	def is_defualt(self):
		"""Tests if the current value of this parameter is equal to the default value
		"""
		return self.default == self.value
	#end is_defualt()
	
	def reset(self):
		"""Resets this parameter to its default value
		"""
		self.value = self.default
	#end reset()
	
# 	def get_value(self):
# 		"""Gets the current value of this parameter, coerced to the type of this parameter
# 		"""
# 		return self._coerce_value(self.value, self.type)
# 	#end get_value()
	
	def get_value_as_type(self, value_type):
		"""Gets the current value of this parameter, coerced to the type specified
		"""
		return self._coerce_value(self.get_value(), value_type)
	#end get_value()
	
	@property
	def type_name(self):
		"""Gets the type of this parameter value as a string representation
		"""
		return self.type.__name__
	#end get_type_as_string()

	def _coerce_value(self, value, value_type):
		new_value = None
		if self.nullable and value is None:
			new_value = None
		elif value_type == bool:
			from distutils.util import strtobool
			new_value = strtobool(str(value))
		elif value_type == int:
			new_value = int(value)
		elif value_type == float:
			new_value = float(value)
		elif value_type == long:
			new_value = long(value)
		elif value_type == complex:
			new_value = complex(value)
		elif value_type == str:
			new_value = str(value)
		elif value_type == list:
			new_value = list(value)
		elif value_type == tuple:
			new_value = tuple(value)
		elif value_type == set:
			new_value = set(value)
		else:
			raise ValueError('Unable to convert %s to type %s for parameter %s!' % (str(value), self.type_name, str(self.name),))
		
		if (self.choices is not None) and (new_value not in self.choices): 
			raise ValueError('The supplied value %s for parameter %s is not present in the allowed choices: %s!' % (str(value), str(self.name), str(self.choices)))
		
		return new_value
	#end _coerce_value()
	
	def __str__(self):
		return "(%s)%s = %s [%s] %s" % (self.type.__name__, str(self.name), str(self.get_value()), str(self.default), str(self.description))
	#end __str__()
	
	def __repr__(self):
		return "ModuleParameter(%s, %s, %s, %s, %s, %s, %s)" % (str(self.name), str(self.type), str(self.default), str(self.value), str(self.description), str(self.nullable), str(self.choices))
	#end __repr__()
#end class ModuleParameter