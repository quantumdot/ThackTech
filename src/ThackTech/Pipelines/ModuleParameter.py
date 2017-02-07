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
		return self.default == self._coerce_value(self.value)
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
	
	def get_type_as_string(self):
		return self.type.__name__
	#end get_type_as_string()

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