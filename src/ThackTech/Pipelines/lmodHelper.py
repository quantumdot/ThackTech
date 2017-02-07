'''
Helper module to support modifying PYTHONPATH via lmod
'''

import logging
import os
import sys
from env_modules_python import module as lmod_module

logger = logging.getLogger(__name__)

def split_path():
	'''
	Split sys.path into the current, system, and pythonpath components
	'''
	# Just the current executable.
	exec_path = sys.path[0:1]

	# Assume that the current PYTHONPATH component of the sys.path is
	# already up to date with the PYTHONPATH environment variable. This
	# helper will FAIL if something else messes with that invariant.
	if 'PYTHONPATH' in os.environ
		python_path = os.environ['PYTHONPATH'].split(":")
	else:
		python_path = sys.path

	# Everything between the executable and the pythonpath.
	# Generally, this is system dependent search paths.
	sys_path = sys.path[1:sys.path.index(python_path[0])]

	# Everything after the pythonpath. Usually, python platform search paths.
	# This is implemented in a weird way to handle the case of the pythonpath
	# including the same entry more than once!
	platform_path = sys.path[sys.path.index(python_path[0])+len(python_path):]

	# Double check our math. Again, this will fail if the invariant is messed with.
	calculated_path = exec_path + sys_path + python_path + platform_path
	if sys.path != calculated_path:
		logger.critical("Calculated sys.path differs from the actual value: expected %s, got %s", sys.path, calculated_path)

	return exec_path, sys_path, python_path, platform_path

def module(subcommand, *args):
	'''
	Wrap the lmod executable in a way that reloads pythonpath
	'''
	logger.debug("Current sys.path: %s", sys.path)

	# Get the soon-to-be-obsolete path information
	exec_path, sys_path, old_python_path, platform_path = split_path()

	# Run the requested lmod command
	logger.debug("Running lmod command: module %s %s", subcommand, " ".join(args)) 
	output = lmod_module(subcommand, *args)

	# TODO: Log a warning if a module has already been imported with the old path?

	# Update the system path
	sys.path = exec_path + sys_path + os.environ['PYTHONPATH'].split(":") + platform_path
	logger.debug("New sys.path: %s", sys.path)

	# Return any output from the lmod exec
	return output